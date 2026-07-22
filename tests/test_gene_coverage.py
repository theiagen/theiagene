import re

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from theiagene import gene_coverage
from theiagene.lib.gene_model import Gene


class MockBam:
    def __init__(self, references, contig_lengths, default_depth=5):
        self.references = tuple(references)
        self._lengths = dict(contig_lengths)
        self._default_depth = default_depth

    def get_reference_length(self, contig):
        return self._lengths[contig]

    def count_coverage(self, contig, start, end, quality_threshold=0):
        span = end - start
        return (
            [self._default_depth] * span,
            [0] * span,
            [0] * span,
            [0] * span,
        )


def _write_mock_gbff(path, contig="contig1", gene="geneA", start=10, end=20):
    record = SeqRecord(Seq("ATCG" * 50), id=contig, name=contig, description="")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(
        SeqFeature(
            FeatureLocation(start, end),
            type="CDS",
            qualifiers={"product": [gene]},
        )
    )
    with open(path, "w") as handle:
        SeqIO.write(record, handle, "genbank")


def test_parse_gbff_extracts_expected_coordinates(tmp_path):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )
    gbff = tmp_path / "mock.gbff"
    _write_mock_gbff(gbff, contig="contig1", gene="geneA", start=10, end=20)

    genes = gene_coverage._coverage_genes(
        gene_coverage.iter_gbff_raw(str(gbff)),
        {"geneA"},
        "product",
        exact_match=True,
        contig_names=set(mock_bam.references),
        require_contig=True,
    )

    assert len(genes) == 1
    assert (genes[0].contig, genes[0].gene_id) == ("contig1", "geneA")
    assert genes[0].cds == [(10, 20)]


def test_bed_and_gbff_coordinates_agree_for_same_gene(tmp_path):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )

    gbff = tmp_path / "mock.gbff"
    _write_mock_gbff(gbff, contig="contig1", gene="geneA", start=10, end=20)

    bed = tmp_path / "mock.bed"
    bed.write_text("contig1\t10\t20\tgeneA\n")

    from_gbff = gene_coverage._coverage_genes(
        gene_coverage.iter_gbff_raw(str(gbff)),
        {"geneA"},
        "product",
        exact_match=True,
        contig_names=set(mock_bam.references),
        require_contig=True,
    )
    from_bed = gene_coverage.parse_bed_genes(
        str(bed),
        {"geneA"},
        exact_match=True,
        contig_names=set(mock_bam.references),
    )

    assert from_gbff[0].cds == from_bed[0].cds == [(10, 20)]


class _Args:
    def __init__(self, **kwargs):
        self.bedfile = None
        self.reference_gbff = None
        self.reference_gff = None
        self.query_genes = None
        self.__dict__.update(kwargs)


def test_input_error_handling_accepts_reference_gff():
    # a GFF reference with query genes is a valid invocation (regression: the
    # GFF option used to be rejected because only GBFF/BED were recognized)
    gene_coverage.input_error_handling(
        _Args(reference_gff="ref.gff", query_genes=["geneA"])
    )


def test_input_error_handling_requires_a_coordinate_source():
    with pytest.raises(FileNotFoundError, match="reference_gff"):
        gene_coverage.input_error_handling(_Args(query_genes=["geneA"]))


def test_input_error_handling_requires_query_genes_without_bed():
    with pytest.raises(ValueError, match="query_genes"):
        gene_coverage.input_error_handling(_Args(reference_gff="ref.gff"))


def test_parse_gff_extracts_expected_coordinates(tmp_path):
    # GFF is 1-based, both-inclusive; geneA at 1-based [11, 20] == 0-based [10, 20)
    gff = tmp_path / "mock.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\t.\tCDS\t11\t20\t.\t+\t0\tID=cds-A;product=geneA\n"
    )
    genes = gene_coverage._coverage_genes(
        gene_coverage.iter_gff_raw(str(gff)),
        {"geneA"},
        "product",
        exact_match=True,
        contig_names={"contig1"},
        require_contig=True,
    )
    assert genes[0].cds == [(10, 20)]


def test_gff_and_bed_coordinates_agree_for_same_gene(tmp_path):
    gff = tmp_path / "mock.gff"
    gff.write_text(
        "##gff-version 3\n"
        "# a comment line that must be skipped\n"
        "contig1\t.\tCDS\t11\t20\t.\t+\t0\tID=cds-A;product=geneA\n"
    )
    bed = tmp_path / "mock.bed"
    bed.write_text("contig1\t10\t20\tgeneA\n")

    from_gff = gene_coverage._coverage_genes(
        gene_coverage.iter_gff_raw(str(gff)),
        {"geneA"}, "product", exact_match=True,
        contig_names={"contig1"}, require_contig=True,
    )
    from_bed = gene_coverage.parse_bed_genes(
        str(bed), {"geneA"}, exact_match=True, contig_names={"contig1"},
    )
    assert from_gff[0].cds == from_bed[0].cds == [(10, 20)]


def test_parse_gff_multi_exon_accumulates_parts(tmp_path):
    # two CDS lines for the same gene accumulate as two coordinate parts
    gff = tmp_path / "multi.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\t.\tCDS\t1\t6\t.\t+\t0\tID=cds-A;product=geneA\n"
        "contig1\t.\tCDS\t11\t16\t.\t+\t0\tID=cds-A;product=geneA\n"
    )
    genes = gene_coverage._coverage_genes(
        gene_coverage.iter_gff_raw(str(gff)),
        {"geneA"}, "product", exact_match=True,
        contig_names={"contig1"}, require_contig=True,
    )
    assert genes[0].cds == [(0, 6), (10, 16)]


def test_parse_gff_raises_when_contig_absent_from_bam(tmp_path):
    gff = tmp_path / "mock.gff"
    gff.write_text(
        "##gff-version 3\n"
        "otherContig\t.\tCDS\t11\t20\t.\t+\t0\tID=cds-A;product=geneA\n"
    )
    with pytest.raises(KeyError, match="otherContig not in BAM"):
        gene_coverage._coverage_genes(
            gene_coverage.iter_gff_raw(str(gff)),
            {"geneA"}, "product", exact_match=True,
            contig_names={"contig1"}, require_contig=True,
        )


def test_quantify_gene_coverage_known_depth_and_breadth():
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 100}, default_depth=5
    )
    genes = [Gene(gene_id="geneA", contig="contig1", parts={"CDS": [(10, 20)]})]

    depth_dict, coverage_dict = gene_coverage.quantify_gene_coverage(
        mock_bam,
        genes,
        min_depth=1,
        min_quality=0,
    )

    assert depth_dict["geneA"] == 5.0
    assert coverage_dict["geneA"] == 100.0


@pytest.mark.parametrize(
    "contig, parts, message_fragment",
    [
        ("contig1", [(10, 10)], "start (10) must be < end (10)"),
        ("contig1", [(5, 15)], "end (15) exceeds contig length (10)"),
        ("missing_contig", [(1, 5)], "not found in BAM references"),
    ],
)
def test_quantify_gene_coverage_edge_guards_raise_clean_value_errors(
    contig, parts, message_fragment
):
    mock_bam = MockBam(
        references=["contig1"], contig_lengths={"contig1": 10}, default_depth=5
    )
    genes = [Gene(gene_id="geneA", contig=contig, parts={"CDS": parts})]

    with pytest.raises(ValueError, match=re.escape(message_fragment)):
        gene_coverage.quantify_gene_coverage(
            mock_bam, genes, min_depth=1, min_quality=0
        )
