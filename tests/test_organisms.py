"""Organism-level regression and quality tests on real TheiaEuk_Illumina_PE data.

Each ``tests/data/<organism>/`` holds a size-minimized subset of one production
sample (see ``tests/data/README.md`` and ``samples.json``): the query genes'
GFF hierarchy, the alignments and VCF records overlapping them, the SnpEff 5.4c
annotation of the extracted variants, and the outputs theiagene produced over
the *full* data (``expected/``).

The subcommands are driven as PHB ``task_gene_coverage``/``task_snpeff`` drive
them, and their outputs are checked three ways:

- regression: the minimized data reproduce the full-data outputs exactly
- provenance: gene coverage matches the production Terra values
- quality: every reported HGVS change is re-derived independently of SnpEff and
  theiagene, from the reference sequence, the GFF CDS model and the organism's
  codon table (table 12 for the CTG-clade *C. auris*), and every reported read
  depth is read back off the VCF sample
"""

import os
import re
import json
import itertools
from urllib.parse import unquote
from collections import defaultdict

import pysam
import pytest

from theiagene import cli
from theiagene import report_variants as rv
from theiagene.lib.parsers import assimilate_gff


DATA = os.path.join(os.path.dirname(__file__), "data")
with open(os.path.join(DATA, "samples.json")) as _handle:
    SAMPLES = json.load(_handle)
ORGANISMS = sorted(SAMPLES)
TERRA_ORGANISMS = sorted(name for name, sample in SAMPLES.items() if "terra" in sample)

# PHB task_gene_coverage/task_snpeff default qualifiers and the medea_magic
# coverage thresholds the production outputs were generated with
FEATURE_QUALIFIER = "product,locus_tag"
COVERAGE_ARGS = ["--min_depth", "10", "--min_mapping_quality", "60", "--min_base_quality", "13"]


def data(organism, name):
    return os.path.join(DATA, organism, name)


def read_tsv(path):
    """Return the rows of a ``#``-headed TSV as dicts"""
    with open(path) as handle:
        header = handle.readline().lstrip("#").rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in handle]


def read_text(path):
    with open(path) as handle:
        return handle.read()


def vcf_sites(path):
    with pysam.VariantFile(path) as handle:
        return [(r.contig, r.pos, r.ref, r.alts) for r in handle]


@pytest.fixture(params=ORGANISMS)
def organism(request):
    return request.param


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """gene_coverage writes its outputs to the working directory"""
    monkeypatch.chdir(tmp_path)
    return tmp_path


def _prepare_database(organism, outdir):
    """Run prepare_snpeff over the organism's GFF; return (genes.gff, config)"""
    template = outdir / "template.config"
    template.write_text("data.dir = ./data/\n")
    config = outdir / "snpEff.config"
    cli.main([
        "prepare_snpeff",
        "--reference_gff", data(organism, "reference.gff"),
        "--template_config", str(template),
        "--data_dir", str(outdir / "data"),
        "--genome_id", "genome",
        "--organism", SAMPLES[organism]["organism"],
        "--output", str(config),
    ])
    return outdir / "data" / "genome" / "genes.gff", config


def _report(organism, gff, output):
    cli.main([
        "report_variants",
        "--vcf", data(organism, "snpeff.vcf"),
        "--reference_gff", str(gff),
        "--query_genes", SAMPLES[organism]["query_genes"],
        "--feature_qualifier", FEATURE_QUALIFIER,
        "--output", str(output),
    ])
    return output


# ---------------------------------------------------------------- gene_coverage


@pytest.fixture
def coverage(organism, workdir):
    """Run gene_coverage as task_gene_coverage does; return the output directory"""
    cli.main([
        "gene_coverage",
        "--bam", data(organism, "reads.bam"),
        "--reference_gff", data(organism, "reference.gff"),
        "--query_genes", SAMPLES[organism]["query_genes"],
        "--feature_qualifier", FEATURE_QUALIFIER,
    ] + COVERAGE_ARGS)
    return workdir


def test_gene_coverage_reproduces_full_data_output(organism, coverage):
    assert read_text(coverage / "COVERAGE_STATS.tsv") == read_text(
        data(organism, "expected/coverage_stats.tsv")
    )


def test_gene_coverage_resolves_and_passes_every_query(organism, coverage):
    """Mirrors task_gene_coverage's QUERY_CHECK: a query is found when its
    quantified length is not 'NA'; production runs hard-fail otherwise"""
    with open(coverage / "LENGTHS_DICT.json") as handle:
        lengths = {k: v for k, v in json.load(handle).items() if k}
    queries = SAMPLES[organism]["query_genes"].split(",")
    assert sorted(lengths) == sorted(queries)
    assert all(length != "NA" and int(length) % 3 == 0 for length in lengths.values())
    for row in read_tsv(coverage / "COVERAGE_STATS.tsv"):
        assert row["reads_mapped_pass"] == "True"
        assert 0 <= float(row["percent_coverage"]) <= 100
        assert float(row["average_depth"]) >= 10


@pytest.mark.parametrize("organism", TERRA_ORGANISMS)
def test_gene_coverage_matches_terra_production(organism, coverage):
    terra = SAMPLES[organism]["terra"]
    rows = {row["query"]: row for row in read_tsv(coverage / "COVERAGE_STATS.tsv")}
    assert sorted(rows) == sorted(terra["gene_coverage_depth_by_gene"])
    for gene, row in rows.items():
        assert float(row["average_depth"]) == pytest.approx(
            float(terra["gene_coverage_depth_by_gene"][gene]), rel=1e-12)
        assert float(row["percent_coverage"]) == float(terra["gene_coverage_breadth_by_gene"][gene])
        assert int(row["reads_mapped"]) == int(terra["gene_coverage_reads_by_gene"][gene])
    # the scalar outputs task_gene_coverage reads back
    assert float((coverage / "MEAN_DEPTH").read_text()) == pytest.approx(
        float(terra["gene_coverage_mean_depth"]), rel=1e-12)
    assert float((coverage / "MEAN_COVERAGE").read_text()) == float(
        terra["gene_coverage_mean_percent_coverage"])
    assert (coverage / "TOTAL_READS").read_text() == terra["gene_coverage_reads_mapped"]


# ------------------------------------------------------------- extract_variants


def test_extract_variants_keeps_only_query_cds_records(organism, workdir):
    """The minimized VCF carries records up to 1 kb beyond the query genes
    (intergenic, intronic, and a neighbouring gene's); extraction must keep
    exactly the records the full-data extraction handed to SnpEff"""
    output = workdir / "extracted.vcf"
    cli.main([
        "extract_variants",
        "--vcf", data(organism, "variants.vcf"),
        "--reference_gff", data(organism, "reference.gff"),
        "--query_genes", SAMPLES[organism]["query_genes"],
        "--feature_qualifier", FEATURE_QUALIFIER,
        "--output", str(output),
    ])
    extracted = vcf_sites(str(output))
    assert extracted == vcf_sites(data(organism, "snpeff.vcf"))
    assert 0 < len(extracted) < len(vcf_sites(data(organism, "variants.vcf")))
    queries = set(SAMPLES[organism]["query_genes"].split(","))
    with pysam.VariantFile(str(output)) as handle:
        assert all(set(record.info["GENE"]) <= queries for record in handle)


# --------------------------------------------------------------- prepare_snpeff


def _attribute_column(path):
    with open(path) as handle:
        return [line.rstrip("\n").split("\t")[8] for line in handle
                if not line.startswith("#") and line.count("\t") == 8]


def test_prepare_snpeff_encodes_pipes_losslessly(organism, tmp_path):
    genes_gff, _ = _prepare_database(organism, tmp_path)
    original = _attribute_column(data(organism, "reference.gff"))
    encoded = _attribute_column(genes_gff)
    assert len(encoded) == len(original)
    assert not any("|" in column for column in encoded)
    # %7C decodes back to the original identifiers, ID/Parent references intact
    assert [column.replace("%7C", "|") for column in encoded] == original


def test_prepare_snpeff_assigns_gff_translation_table(organism, tmp_path):
    """C. auris is a CTG-clade yeast (CUG -> Ser); its GFF declares
    transl_table=12 per CDS, which must reach SnpEff per contig. The others
    declare none and fall back to SnpEff's Standard table"""
    _, config = _prepare_database(organism, tmp_path)
    codon_lines = [line for line in config.read_text().splitlines() if "codonTable" in line]
    cds_contigs = {
        line.split("\t")[0] for line in open(data(organism, "reference.gff"))
        if "\tCDS\t" in line
    }
    if SAMPLES[organism]["translation_table"] == 12:
        assert sorted(codon_lines) == sorted(
            f"genome.{contig}.codonTable : Alternative_Yeast_Nuclear" for contig in cds_contigs
        )
    else:
        assert codon_lines == []


# -------------------------------------------------------------- report_variants


@pytest.mark.parametrize("gff_source", ["reference", "database"])
def test_report_variants_reproduces_full_data_output(organism, gff_source, tmp_path):
    """Production reports against prepare_snpeff's percent-encoded genes.gff;
    the original GFF must resolve the same records"""
    if gff_source == "database":
        gff, _ = _prepare_database(organism, tmp_path)
    else:
        gff = data(organism, "reference.gff")
    report = _report(organism, gff, tmp_path / "report.tsv")
    assert read_text(report) == read_text(data(organism, "expected/variant_report.tsv"))


def test_report_variants_rejects_production_annotation_with_unencoded_pipes():
    """The Terra annotation of C. auris predates prepare_snpeff's pipe encoding:
    NCBI's 'rna-gnl|WGS:PGLS|mrna_...' transcript IDs shifted every later ANN
    subfield, so Terra reported 'protein_coding'/'1/1' as the HGVS changes.
    Such an annotation must raise rather than be misreported"""
    organism = "candidozyma_auris"
    features = assimilate_gff(data(organism, "reference.gff"))
    with pytest.raises(ValueError, match="prepare_snpeff"):
        rv.report_variants(
            data(organism, "snpeff_unencoded.vcf"), features,
            {"upstream_gene_variant", "downstream_gene_variant"}, "CDS",
            FEATURE_QUALIFIER.split(","), SAMPLES[organism]["query_genes"].split(","),
        )


def test_report_variants_recovers_hgvs_shifted_by_unencoded_pipes():
    """The same two C. auris variants, re-annotated from the encoded GFF"""
    rows = read_tsv(data("candidozyma_auris", "expected/variant_report.tsv"))
    assert [(row["NT"], row["AA"]) for row in rows] == [
        ("1113dupG", "R372fs"), ("807_809delGCA", "Q270del"),
    ]
    assert all(
        row["HGVSc"].startswith('"rna-gnl|WGS:PGLS|mrna_CJI97_000402T0":c.') for row in rows
    )


def test_report_variants_cyp51A_reports_known_polymorphism_sites():
    """A. fumigatus cyp51A (AFUA_4G06890) varies at the well-characterised
    F46Y/M172V/N248T/D255E/E427K sites; the Af293 reference carries the
    Y/V/T/E/K residues there, so this isolate reads as their reverse"""
    rows = read_tsv(data("aspergillus_fumigatus", "expected/variant_report.tsv"))
    cyp51a = {row["AA"] for row in rows if row["GENE"] == "AFUA_4G06890"}
    assert {"Y46F", "V172M", "T248N", "E255D", "K427E"} <= cyp51a


# independent re-derivation of the reported HGVS changes

_BASES = "TCAG"
_STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def codon_table(number):
    table = {"".join(codon): aa for codon, aa in zip(itertools.product(_BASES, repeat=3), _STANDARD)}
    if number == 12:
        table["CTG"] = "S"
    elif number != 1:
        raise ValueError(f"no codon table {number} in the test helper")
    return table


def load_reference(path):
    """Return ``fetch(contig, start, end)`` (1-based, inclusive) over the
    coordinate-labelled ``<contig>:<start>-<end>`` span records"""
    spans = []
    with pysam.FastxFile(path) as handle:
        for record in handle:
            contig, coords = record.name.rsplit(":", 1)
            lo, hi = map(int, coords.split("-"))
            spans.append((contig, lo, hi, record.sequence.upper()))

    def fetch(contig, start, end):
        for name, lo, hi, seq in spans:
            if name == contig and lo <= start and end <= hi:
                return seq[start - lo:end - lo + 1]
        raise KeyError(f"{contig}:{start}-{end} outside the reference spans")

    return fetch


def load_cds(path):
    """Return ``{transcript ID: [(contig, start, end, strand, phase), ...]}``
    read directly off the GFF's CDS records"""
    cds = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = dict(pair.split("=", 1) for pair in fields[8].split(";") if pair)
            for parent in attrs["Parent"].split(","):
                cds[unquote(parent)].append(
                    (fields[0], int(fields[3]), int(fields[4]), fields[6], fields[7])
                )
    return cds


def coding_sequence(segments, fetch):
    """Splice a transcript's CDS in transcription order"""
    strand = segments[0][3]
    ordered = sorted(segments, key=lambda s: s[1], reverse=strand == "-")
    if ordered[0][4] not in ("0", "."):
        return None  # 5'-partial CDS: codon frame is not anchored at c.1
    seq = "".join(fetch(c, s, e) for c, s, e, _, _ in sorted(segments, key=lambda s: s[1]))
    return seq.translate(_COMPLEMENT)[::-1] if strand == "-" else seq


def split_prefixed(hgvs):
    """Split ``identifier:hgvs`` / ``"identifier:with:colons":hgvs``"""
    if hgvs.startswith('"'):
        identifier, change = hgvs[1:].split('":', 1)
    else:
        identifier, change = hgvs.split(":", 1)
    return identifier, change


_CDNA = re.compile(
    r"^c\.(\d+)(?:_(\d+))?"
    r"(?:([ACGT])>([ACGT])|del([ACGT]+)(?:ins([ACGT]+))?|dup([ACGT]+))$"
)
_PROTEIN = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)(?:([A-Z][a-z]{2}|\*)$)?")
# SnpEff 5.4c writes a synonymous multi-codon delins as a bare 'p.<N>'
_BARE_PROTEIN = re.compile(r"^p\.(\d+)$")


def test_report_variants_hgvs_agrees_with_reference_sequence(organism):
    """Re-derive each reported change from the reference: the HGVS.c reference
    bases must be the transcript's CDS bases at those positions and the HGVS.p
    reference residue the translation of that codon under the organism's codon
    table. For a substitution (SNV, or a same-length delins/MNP) the affected
    codons are re-translated with the alternate bases: a missense/stop change
    must alter exactly the reported residue, a synonymous one none. Every record
    must be checkable (extraction keeps only CDS variants, so none is intronic)"""
    fetch = load_reference(data(organism, "reference.fasta"))
    cds = load_cds(data(organism, "reference.gff"))
    table = codon_table(SAMPLES[organism]["translation_table"])
    rows = read_tsv(data(organism, "expected/variant_report.tsv"))
    assert rows

    for row in rows:
        transcript, hgvsc = split_prefixed(row["HGVSc"])
        seq = coding_sequence(cds[transcript], fetch)
        assert seq is not None, f"{transcript} CDS is 5'-partial"
        assert len(seq) % 3 == 0 and table[seq[-3:]] == "*"
        cdna = _CDNA.match(hgvsc)
        assert cdna, f"unrecognised HGVS.c {hgvsc}"
        first, last, sub_ref, sub_alt, deleted, inserted, duplicated = cdna.groups()
        start, end = int(first), int(last or first)
        ref_bases = sub_ref or deleted or duplicated
        assert seq[start - 1:end] == ref_bases, f"{row['HGVSc']} vs CDS {seq[start - 1:end]}"

        # residue-level change over the affected codons, for same-length changes
        alt_bases = sub_alt or (inserted if inserted and len(inserted) == len(deleted) else None)
        changed = None
        if alt_bases:
            mutated = seq[:start - 1] + alt_bases + seq[end:]
            codons = range((start - 1) // 3, (end - 1) // 3 + 1)
            changed = [
                (i + 1, table[seq[3 * i:3 * i + 3]], table[mutated[3 * i:3 * i + 3]])
                for i in codons
                if table[seq[3 * i:3 * i + 3]] != table[mutated[3 * i:3 * i + 3]]
            ]

        if row["HGVSp"] == "NA":
            continue
        _, hgvsp = split_prefixed(row["HGVSp"])
        bare = _BARE_PROTEIN.match(hgvsp)
        if bare:
            # undocumented SnpEff notation: the change must be synonymous across
            # every codon, and N is one past the last affected codon (pinned so
            # a change in SnpEff's or theiagene's handling surfaces here)
            consequence = row["REPORT"].split(" (", 1)[1].split(" ", 1)[0]
            assert "synonymous_variant" in consequence.split("&")
            assert changed == [], row["HGVSc"]
            assert int(bare.group(1)) == codons[-1] + 2
            continue
        protein = _PROTEIN.match(hgvsp)
        assert protein, f"unrecognised HGVS.p {hgvsp}"
        ref_aa, aa_pos, alt_aa = protein.groups()
        aa_pos = int(aa_pos)
        assert table[seq[3 * aa_pos - 3:3 * aa_pos]] == _THREE_TO_ONE[ref_aa], row["HGVSp"]
        if changed is not None:
            assert aa_pos - 1 in codons
            if alt_aa == ref_aa:
                assert changed == [], row["HGVSp"]
            else:
                alt_one = "*" if alt_aa == "*" else _THREE_TO_ONE[alt_aa]
                assert changed == [(aa_pos, _THREE_TO_ONE[ref_aa], alt_one)], row["HGVSp"]


def test_report_variants_columns_are_consistent(organism):
    """GENE is the query term asked for, NT/AA abbreviate HGVSc/HGVSp, the
    REPORT line leads with them, and its read depths are the VCF sample's
    (AD for GATK; RO/AO for the FreeBayes-called C. neoformans)"""
    queries = SAMPLES[organism]["query_genes"].split(",")
    depths = {}
    with pysam.VariantFile(data(organism, "snpeff.vcf")) as handle:
        for record in handle:
            sample = next(iter(record.samples.values()))
            if "AD" in sample.keys() and sample["AD"][0] is not None:
                ref_depth, alt_depths = sample["AD"][0], sample["AD"][1:]
            else:
                ref_depth, alt_depths = sample["RO"], sample["AO"]
            for alt, depth in zip(record.alts, alt_depths):
                for raw in record.info["ANN"]:
                    fields = raw.split("|")
                    if fields[0] == alt and fields[9]:
                        depths[(unquote(fields[6]), fields[9])] = f"{record.ref}:{ref_depth} {alt}:{depth}"

    rows = read_tsv(data(organism, "expected/variant_report.tsv"))
    for row in rows:
        transcript, hgvsc = split_prefixed(row["HGVSc"])
        assert row["GENE"] in queries
        assert row["NT"] == hgvsc[2:]
        if row["HGVSp"] != "NA":
            hgvsp = split_prefixed(row["HGVSp"])[1]
            abbreviated = re.sub("|".join(_THREE_TO_ONE), lambda m: _THREE_TO_ONE[m.group()], hgvsp[2:])
            assert row["AA"] == abbreviated
            assert f" {hgvsc} {hgvsp}; " in row["REPORT"]
        assert row["REPORT"].startswith(f'{row["GENE"]}: "')
        assert row["REPORT"].endswith(f"; {depths[(transcript, hgvsc)]})")
