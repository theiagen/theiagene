#!/usr/bin/env python3
"""Rebuild the size-minimized organism test data under ``tests/data/``.

Not run by the test suite. This records how the committed fixtures were derived
from production TheiaEuk_Illumina_PE outputs so they can be regenerated:

1. download each sample's BAM, VCF and reference FASTA/GFF (``samples.json``
   ``sources``; ``gsutil`` + Terra workspace read access required)
2. run theiagene (this checkout's ``src/``) over the *full* data in the same
   order as PHB ``task_gene_coverage``/``task_snpeff``: gene_coverage;
   prepare_snpeff -> extract_variants -> SnpEff build/ann -> report_variants.
   These full-data outputs become ``expected/``
3. minimize each input to the query genes: the GFF features overlapping the
   query-gene spans, the alignments overlapping those spans (aux tags dropped,
   read names shortened), the VCF records within ``VCF_FLANK`` of them, and the
   span sequences as a coordinate-labelled FASTA (``<contig>:<start>-<end>``)
4. re-run gene_coverage/report_variants on the minimized data and fail unless
   they reproduce the full-data outputs exactly

Requirements: pysam, ``samtools`` >= 1.10 and docker with a SnpEff image (only
its ``snpeff`` jar wrapper and bundled ``snpEff.config`` are used).

    python3 tests/data/build_test_data.py --work /path/to/scratch \\
        --samtools /path/to/samtools \\
        --snpeff_image us-docker.pkg.dev/general-theiagen/theiagen/snpeff:5.4c
"""

import os
import re
import sys
import gzip
import json
import shutil
import argparse
import subprocess

import pysam

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(DATA_DIR))
SRC = os.path.join(REPO, "src")
sys.path.insert(0, SRC)

from theiagene.lib.parsers import assimilate_gff  # noqa: E402
from theiagene.lib.query import gff_query_ranges, ordered_query_genes, split_qualifiers  # noqa: E402

# VCF records this far beyond a query-gene span are kept, so extraction has
# intergenic/intronic records to exclude
VCF_FLANK = 1000
# PHB medea_magic -> task_gene_coverage parameters
COVERAGE_ARGS = ["--min_depth", "10", "--min_mapping_quality", "60", "--min_base_quality", "13"]
# PHB task_gene_coverage / task_snpeff default
FEATURE_QUALIFIER = "product,locus_tag"
# header lines that only carry execution paths/commands
_DROP_HEADER = re.compile(r"^##(GATKCommandLine|bcftools_|commandline|reference=|fileDate)")


def run(cmd, **kwargs):
    print("+", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, **kwargs)


def theiagene(args, cwd):
    env = dict(os.environ, PYTHONPATH=SRC)
    run([sys.executable, "-m", "theiagene"] + args, cwd=cwd, env=env)


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def download(sources: dict, raw_dir: str) -> dict:
    os.makedirs(raw_dir, exist_ok=True)
    local = {}
    for key, uri in sources.items():
        path = os.path.join(raw_dir, os.path.basename(uri))
        if not os.path.exists(path):
            run(["gsutil", "-q", "cp", uri, path])
        local[key] = path
    return local


def fasta_lengths(fasta: str) -> dict:
    pysam.faidx(fasta)
    with pysam.FastaFile(fasta) as handle:
        return dict(zip(handle.references, handle.lengths))


def normalize_vcf(vcf: str, out: str, rename: dict, strip_ann: bool, lengths: dict) -> str:
    """Rename contigs and drop a prior ``ANN`` annotation (e.g. Snippy's) so the
    VCF matches the reference and can be re-annotated."""
    contigs_declared = False
    with open_text(vcf) as src, open(out, "w") as dst:
        for line in src:
            if line.startswith("##"):
                if strip_ann and line.startswith(("##INFO=<ID=ANN,", "##SnpEff")):
                    continue
                if line.startswith("##contig=<ID="):
                    contigs_declared = True
                    for old, new in rename.items():
                        line = line.replace(f"<ID={old},", f"<ID={new},")
                dst.write(line)
                continue
            if line.startswith("#CHROM"):
                if not contigs_declared:
                    dst.writelines(f"##contig=<ID={c},length={n}>\n" for c, n in lengths.items())
                dst.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            fields[0] = rename.get(fields[0], fields[0])
            if strip_ann:
                info = [kv for kv in fields[7].split(";") if not kv.startswith("ANN=")]
                fields[7] = ";".join(info) or "."
            dst.write("\t".join(fields) + "\n")
    return out


def query_spans(gff: str, query_genes: str) -> dict:
    """Return ``{contig: [(start, end), ...]}`` (1-based, inclusive) covering
    every query gene's full GFF hierarchy, expanded to whole overlapping genes."""
    features = assimilate_gff(gff)
    ranges = gff_query_ranges(
        features, ordered_query_genes(query_genes), "CDS",
        split_qualifiers(FEATURE_QUALIFIER), False,
    )
    by_label = {}
    for contig, segments in ranges.items():
        for start, end, label in segments:
            lo, hi = by_label.get((contig, label), (start + 1, end))
            by_label[(contig, label)] = (min(lo, start + 1), max(hi, end))
    spans = {}
    for (contig, _), span in by_label.items():
        spans.setdefault(contig, []).append(span)

    # expand each span to the full extent of every gene overlapping it until stable
    genes = []
    with open_text(gff) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            fields = line.split("\t")
            if len(fields) == 9 and fields[2] == "gene":
                genes.append((fields[0], int(fields[3]), int(fields[4])))
    changed = True
    while changed:
        changed = False
        for contig, intervals in spans.items():
            merged = _merge(intervals)
            for i, (lo, hi) in enumerate(merged):
                for gcontig, gstart, gend in genes:
                    if gcontig == contig and gstart <= hi and gend >= lo:
                        if gstart < lo or gend > hi:
                            lo, hi = min(lo, gstart), max(hi, gend)
                            changed = True
                merged[i] = (lo, hi)
            spans[contig] = _merge(merged)
    return spans


def _merge(intervals):
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def _overlaps(spans, contig, start, end, flank=0):
    return any(lo - flank <= end and hi + flank >= start for lo, hi in spans.get(contig, ()))


def minimize_gff(gff: str, spans: dict, lengths: dict, out: str) -> None:
    with open_text(gff) as src, open(out, "w") as dst:
        dst.write("##gff-version 3\n")
        for contig in spans:
            dst.write(f"##sequence-region {contig} 1 {lengths[contig]}\n")
        for line in src:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9 or fields[0] not in spans:
                continue
            if fields[2] == "region" or _overlaps(spans, fields[0], int(fields[3]), int(fields[4])):
                dst.write(line)


def minimize_fasta(fasta: str, spans: dict, out: str) -> None:
    with pysam.FastaFile(fasta) as handle, open(out, "w") as dst:
        for contig, intervals in spans.items():
            for lo, hi in intervals:
                seq = handle.fetch(contig, lo - 1, hi)
                dst.write(f">{contig}:{lo}-{hi}\n")
                dst.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))


def minimize_bam(bam: str, spans: dict, out: str) -> None:
    """Write the alignments overlapping ``spans`` with aux tags dropped and read
    names shortened (consistently across mates); the header keeps only @HD/@SQ."""
    names = {}
    seen = set()
    unsorted = out + ".unsorted.bam"
    with pysam.AlignmentFile(bam) as src:
        header = src.header.to_dict()
        header = {"HD": header.get("HD", {"VN": "1.6"}), "SQ": header["SQ"]}
        with pysam.AlignmentFile(unsorted, "wb", header=header) as dst:
            for contig, intervals in spans.items():
                for lo, hi in intervals:
                    for read in src.fetch(contig, lo - 1, hi):
                        key = (read.query_name, read.flag, read.reference_start)
                        if key in seen:
                            continue
                        seen.add(key)
                        name = names.setdefault(read.query_name, f"r{len(names) + 1}")
                        new = pysam.AlignedSegment.from_dict(read.to_dict(), dst.header)
                        new.query_name = name
                        new.set_tags([])
                        dst.write(new)
    pysam.sort("-o", out, unsorted)
    os.remove(unsorted)
    pysam.index(out)


def minimize_vcf(vcf: str, spans: dict, out: str, flank: int) -> int:
    kept = 0
    with open_text(vcf) as src, open(out, "w") as dst:
        for line in src:
            if line.startswith("##"):
                if _DROP_HEADER.match(line):
                    continue
                match = re.match(r"##contig=<ID=([^,>]+)", line)
                if match and match.group(1) not in spans:
                    continue
                dst.write(line)
            elif line.startswith("#"):
                dst.write(line)
            else:
                fields = line.split("\t", 5)
                pos = int(fields[1])
                if _overlaps(spans, fields[0], pos, pos + len(fields[3]) - 1, flank):
                    dst.write(line)
                    kept += 1
    return kept


def sanitize_snpeff_vcf(vcf: str, out: str) -> None:
    """Copy a SnpEff VCF minus execution-path header lines."""
    with open(vcf) as src, open(out, "w") as dst:
        dst.writelines(line for line in src if not _DROP_HEADER.match(line))


def run_gene_coverage(bam: str, gff: str, query_genes: str, outdir: str) -> str:
    os.makedirs(outdir, exist_ok=True)
    theiagene(["gene_coverage", "--bam", os.path.abspath(bam), "--reference_gff",
               os.path.abspath(gff), "--query_genes", query_genes,
               "--feature_qualifier", FEATURE_QUALIFIER] + COVERAGE_ARGS, cwd=outdir)
    return os.path.join(outdir, "COVERAGE_STATS.tsv")


def run_snpeff(cfg: dict, raw: dict, vcf: str, work: str, image: str) -> tuple:
    """Mirror PHB task_snpeff over the full data, returning the annotated VCF
    and the report TSV."""
    genome_id = re.sub(r"[^A-Za-z0-9._-]", "_", cfg["sample"])
    snp_dir = os.path.join(work, "snpeff")
    os.makedirs(snp_dir, exist_ok=True)
    template = os.path.join(snp_dir, "template.config")
    with open(template, "w") as handle:
        subprocess.run(["docker", "run", "--rm", image, "cat", "/snpEff/snpEff.config"],
                       check=True, stdout=handle)
    theiagene(["prepare_snpeff", "--reference_gff", raw["gff"], "--template_config",
               template, "--data_dir", "data", "--genome_id", genome_id, "--organism",
               cfg["organism"], "--output", "snpEff.config"], cwd=snp_dir)
    shutil.copy(raw["fasta"], os.path.join(snp_dir, "data", genome_id, "sequences.fa"))
    db_gff = os.path.join(snp_dir, "data", genome_id, "genes.gff")
    extracted = f"{genome_id}_extracted.vcf"
    theiagene(["extract_variants", "--vcf", vcf, "--reference_gff", db_gff,
               "--query_genes", cfg["query_genes"], "--feature_qualifier",
               FEATURE_QUALIFIER, "--output", extracted], cwd=snp_dir)
    uid = f"{os.getuid()}:{os.getgid()}"
    run(["docker", "run", "--rm", "-u", uid, "-v", f"{snp_dir}:/work", "-w", "/work",
         "-e", "JAVA_TOOL_OPTIONS=-Xmx8G", image, "bash", "-c",
         # data.dir points at the host path, so re-point it inside the container
         f"sed -i 's#^data.dir = .*#data.dir = /work/data/#' snpEff.config && "
         f"snpeff build -c snpEff.config -gff3 -noCheckCds -noCheckProtein -noLog {genome_id} && "
         f"snpeff ann -c snpEff.config -nodownload -noLog -ud 0 -stats summary.html "
         f"{genome_id} {extracted} > annotated.vcf"])
    annotated = os.path.join(snp_dir, "annotated.vcf")
    report = os.path.join(snp_dir, "variant_report.tsv")
    theiagene(["report_variants", "--vcf", annotated, "--reference_gff", db_gff,
               "--query_genes", cfg["query_genes"], "--feature_qualifier",
               FEATURE_QUALIFIER, "--output", report], cwd=snp_dir)
    return annotated, report


def build(name: str, cfg: dict, args) -> None:
    work = os.path.join(args.work, name)
    out = os.path.join(DATA_DIR, name)
    os.makedirs(os.path.join(out, "expected"), exist_ok=True)
    raw = download(cfg["sources"], os.path.join(work, "raw"))
    lengths = fasta_lengths(raw["fasta"])

    # reconcile contig names (e.g. Snippy's unversioned 'NC_006670') with the reference
    rename = cfg.get("contig_rename", {})
    vcf = raw["vcf"]
    bam = raw["bam"]
    if rename or cfg.get("strip_ann"):
        vcf = normalize_vcf(raw["vcf"], os.path.join(work, "normalized.vcf"),
                            rename, cfg.get("strip_ann", False), lengths)
    if rename:
        header = os.path.join(work, "header.sam")
        with open(header, "w") as handle:
            text = subprocess.run([args.samtools, "view", "-H", raw["bam"]], check=True,
                                  capture_output=True, text=True).stdout
            for old, new in rename.items():
                text = text.replace(f"SN:{old}\t", f"SN:{new}\t")
            handle.write(text)
        bam = os.path.join(work, "renamed.bam")
        with open(bam, "wb") as handle:
            subprocess.run([args.samtools, "reheader", header, raw["bam"]], check=True, stdout=handle)
        run([args.samtools, "index", bam])

    # full-data outputs
    full_cov = run_gene_coverage(bam, raw["gff"], cfg["query_genes"], os.path.join(work, "full_coverage"))
    annotated, report = run_snpeff(cfg, raw, os.path.abspath(vcf), work, args.snpeff_image)

    # minimized inputs
    spans = query_spans(raw["gff"], cfg["query_genes"])
    minimize_gff(raw["gff"], spans, lengths, os.path.join(out, "reference.gff"))
    minimize_fasta(raw["fasta"], spans, os.path.join(out, "reference.fasta"))
    minimize_bam(bam, spans, os.path.join(out, "reads.bam"))
    kept = minimize_vcf(vcf, spans, os.path.join(out, "variants.vcf"), VCF_FLANK)
    sanitize_snpeff_vcf(annotated, os.path.join(out, "snpeff.vcf"))
    shutil.copy(full_cov, os.path.join(out, "expected", "coverage_stats.tsv"))
    shutil.copy(report, os.path.join(out, "expected", "variant_report.tsv"))
    if "terra_snpeff_vcf" in raw:
        # the production annotation predating prepare_snpeff's pipe encoding
        sanitize_snpeff_vcf(raw["terra_snpeff_vcf"], os.path.join(out, "snpeff_unencoded.vcf"))
    print(f"{name}: {sum(map(len, spans.values()))} span(s), {kept} VCF record(s)", file=sys.stderr)

    # the minimized data must reproduce the full-data outputs
    mini_cov = run_gene_coverage(os.path.join(out, "reads.bam"), os.path.join(out, "reference.gff"),
                                 cfg["query_genes"], os.path.join(work, "mini_coverage"))
    _assert_same(full_cov, mini_cov)
    mini_report = os.path.join(work, "mini_variant_report.tsv")
    theiagene(["report_variants", "--vcf", os.path.join(out, "snpeff.vcf"), "--reference_gff",
               os.path.join(out, "reference.gff"), "--query_genes", cfg["query_genes"],
               "--feature_qualifier", FEATURE_QUALIFIER, "--output", mini_report], cwd=work)
    _assert_same(report, mini_report)


def _assert_same(expected: str, observed: str) -> None:
    with open(expected) as a, open(observed) as b:
        if a.read() != b.read():
            raise SystemExit(f"minimized output {observed} differs from full-data {expected}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--work", required=True, help="scratch directory for downloads/full-data runs")
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--snpeff_image", default="us-docker.pkg.dev/general-theiagen/theiagen/snpeff:5.4c")
    parser.add_argument("--organisms", nargs="*", help="DEFAULT: all in samples.json")
    args = parser.parse_args(argv)
    with open(os.path.join(DATA_DIR, "samples.json")) as handle:
        samples = json.load(handle)
    for name in args.organisms or samples:
        build(name, samples[name], args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
