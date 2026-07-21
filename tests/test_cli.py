"""End-to-end tests for the unified ``theiagene`` command-line entrypoint."""

import pytest
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from theiagene import cli, gene_coverage, variant_annotation


ALPHA_CODING = "ATGTATCCCAAAGGGTTTCATTGA"  # M Y P K G F H *
ALPHA_START = 9


def _write_alpha_gbff(path):
    """Minimal single forward-strand gene GBFF (mirrors the annotation tests)."""
    rec = SeqRecord(
        Seq("A" * ALPHA_START + ALPHA_CODING + "A" * 7),
        id="chr1",
        name="chr1",
        description="",
    )
    rec.annotations["molecule_type"] = "DNA"
    rec.features.append(
        SeqFeature(
            FeatureLocation(ALPHA_START, ALPHA_START + len(ALPHA_CODING), strand=1),
            type="CDS",
            qualifiers={"product": ["test gene alpha"], "transl_table": ["1"]},
        )
    )
    with open(path, "w") as handle:
        SeqIO.write([rec], handle, "genbank")


def _write_missense_vcf(path):
    lines = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=chr1,length=40>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t14\t.\tA\tT\t.\t.\t.",  # c.5 A>T, codon 2 TAT->TTT (Tyr2Phe)
    ]
    path.write_text("\n".join(lines) + "\n")


def _write_bam(path, contig="contig1", contig_len=100, read_start=10, read_len=50):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": contig_len, "SN": contig}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        seg = pysam.AlignedSegment()
        seg.query_name = "read1"
        seg.query_sequence = "A" * read_len
        seg.flag = 0
        seg.reference_id = 0
        seg.reference_start = read_start
        seg.mapping_quality = 60
        seg.cigar = [(0, read_len)]  # all match
        seg.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        out.write(seg)
    pysam.index(str(path))


# --------------------------------------------------------------------------- #
# wiring
# --------------------------------------------------------------------------- #

def test_build_parser_registers_both_subcommands():
    parser = cli.build_parser()

    gc = parser.parse_args(["gene_coverage", "--bam", "x.bam", "--bedfile", "y.bed"])
    assert gc._handler is gene_coverage.run_cli
    assert gc.bam == "x.bam"

    va = parser.parse_args(
        ["variant_annotation", "--vcf", "v.vcf", "--reference_gbff", "r.gbff",
         "--query_genes", "geneA", "geneB"]
    )
    assert va._handler is variant_annotation.run_cli
    assert va.query_genes == ["geneA", "geneB"]


def test_version_flag_exits_cleanly(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0
    assert "theiagene" in capsys.readouterr().out


def test_missing_subcommand_errors():
    with pytest.raises(SystemExit) as exc:
        cli.main([])
    assert exc.value.code != 0


# --------------------------------------------------------------------------- #
# variant_annotation subcommand, end to end
# --------------------------------------------------------------------------- #

def test_variant_annotation_subcommand_end_to_end(tmp_path, capsys):
    gbff = tmp_path / "ref.gbff"
    vcf = tmp_path / "v.vcf"
    out = tmp_path / "annotations.txt"
    _write_alpha_gbff(gbff)
    _write_missense_vcf(vcf)

    rc = cli.main([
        "variant_annotation",
        "--vcf", str(vcf),
        "--reference_gbff", str(gbff),
        "--query_genes", "test gene alpha",
        "--exact_match",
        "--output", str(out),
    ])

    assert rc == 0
    report = out.read_text()
    assert "missense_variant c.5A>T p.Tyr2Phe" in report
    assert "missense_variant c.5A>T p.Tyr2Phe" in capsys.readouterr().out


def test_variant_annotation_requires_a_reference(tmp_path):
    vcf = tmp_path / "v.vcf"
    _write_missense_vcf(vcf)
    with pytest.raises(ValueError, match="reference_gbff"):
        cli.main([
            "variant_annotation",
            "--vcf", str(vcf),
            "--query_genes", "test gene alpha",
        ])


# --------------------------------------------------------------------------- #
# gene_coverage subcommand, end to end
# --------------------------------------------------------------------------- #

def test_gene_coverage_subcommand_end_to_end(tmp_path, monkeypatch):
    bam = tmp_path / "reads.bam"
    bed = tmp_path / "regions.bed"
    _write_bam(bam, contig="contig1", contig_len=100, read_start=10, read_len=50)
    bed.write_text("contig1\t10\t20\tgeneA\n")

    # run_cli writes its outputs into the current working directory
    monkeypatch.chdir(tmp_path)
    rc = cli.main(["gene_coverage", "--bam", str(bam), "--bedfile", str(bed)])

    assert rc == 0
    tsv = (tmp_path / "COVERAGE_STATS.tsv").read_text()
    assert "geneA" in tsv
    # the single 50bp read fully spans geneA's [10, 20) at depth 1
    assert "geneA\t1.0\t100.0" in tsv
    assert (tmp_path / "DEPTH_DICT.json").exists()
    assert (tmp_path / "COVERAGE_DICT.json").exists()
