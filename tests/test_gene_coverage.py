"""Unit tests for theiagene.gene_coverage (coverage quantification + helpers).

Coverage is quantified over a ``{<CONTIG>: [(START, END, LABEL), ...]}`` map (the
shape :func:`theiagene.lib.query.gff_query_ranges`/``bed_query_ranges`` produce);
every range sharing a ``LABEL`` is pooled onto one output row."""

import re

import pytest

from theiagene import gene_coverage


# --------------------------------------------------------------------------- #
# test doubles
# --------------------------------------------------------------------------- #

class MockBam:
    """A minimal stand-in for pysam.AlignmentFile.

    ``depth_fn(contig, pos) -> int`` supplies per-base depth (default: a flat
    ``default_depth``); all of it is reported in the first (A) channel so the
    four-channel sum used by quantify_gene_coverage equals that depth."""

    def __init__(self, lengths, depth_fn=None, default_depth=5):
        self._lengths = dict(lengths)
        self.references = tuple(lengths)
        self._depth_fn = depth_fn
        self._default = default_depth
        self.seen_quality = []

    def get_reference_length(self, contig):
        return self._lengths[contig]

    def count_coverage(self, contig, start, end, quality_threshold=0):
        self.seen_quality.append(quality_threshold)
        a = [
            self._depth_fn(contig, pos) if self._depth_fn else self._default
            for pos in range(start, end)
        ]
        zeros = [0] * (end - start)
        return (a, list(zeros), list(zeros), list(zeros))


class _Args:
    """A stand-in argparse.Namespace for input_error_handling."""

    def __init__(self, **kwargs):
        self.bedfile = None
        self.reference_gff = None
        self.query_genes = None
        self.__dict__.update(kwargs)


# --------------------------------------------------------------------------- #
# input_error_handling
# --------------------------------------------------------------------------- #

def test_input_error_handling_accepts_gff_with_query_genes():
    gene_coverage.input_error_handling(
        _Args(reference_gff="ref.gff", query_genes=["geneA"])
    )


def test_input_error_handling_accepts_bedfile_alone():
    # a BED supplies both a coordinate source and the query names
    gene_coverage.input_error_handling(_Args(bedfile="regions.bed"))


def test_input_error_handling_requires_a_coordinate_source():
    with pytest.raises(FileNotFoundError, match="reference_gff.*bedfile|bedfile"):
        gene_coverage.input_error_handling(_Args(query_genes=["geneA"]))


def test_input_error_handling_requires_query_genes_without_bed():
    with pytest.raises(ValueError, match="query_genes"):
        gene_coverage.input_error_handling(_Args(reference_gff="ref.gff"))


# --------------------------------------------------------------------------- #
# make_tsv
# --------------------------------------------------------------------------- #

def test_make_tsv_renders_header_and_rows():
    tsv = gene_coverage.make_tsv(
        {"geneA": 5.0, "geneB": 2.0},
        {"geneA": 100.0, "geneB": 50.0},
        ambiguous_contig=False,
    )
    lines = tsv.splitlines()
    assert lines[0] == "#query\taverage_depth\tpercent_coverage"
    assert lines[1] == "geneA\t5.0\t100.0"
    assert lines[2] == "geneB\t2.0\t50.0"


def test_make_tsv_warns_in_header_when_contig_is_ambiguous():
    tsv = gene_coverage.make_tsv({"geneA": 1.0}, {"geneA": 100.0}, ambiguous_contig=True)
    assert tsv.splitlines()[0].startswith("#query (WARNING")


# --------------------------------------------------------------------------- #
# quantify_gene_coverage
# --------------------------------------------------------------------------- #

def test_quantify_known_depth_and_full_breadth():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    assert depth["geneA"] == 5.0
    assert coverage["geneA"] == 100.0


def test_quantify_pools_segments_sharing_a_label():
    bam = MockBam({"contig1": 100}, default_depth=4)
    # two CDS segments of one query gene share a label and are pooled
    ranges = {"contig1": [(0, 6, "geneA"), (10, 16, "geneA")]}

    depth, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    # 12 evaluated bases, all at depth 4
    assert depth["geneA"] == 4.0
    assert coverage["geneA"] == 100.0


def test_quantify_partial_breadth_from_min_depth_threshold():
    # first half of the region at depth 10, second half uncovered
    def depth_fn(contig, pos):
        return 10 if pos < 15 else 0

    bam = MockBam({"contig1": 100}, depth_fn=depth_fn)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    assert depth["geneA"] == 5.0     # mean of five 10s and five 0s
    assert coverage["geneA"] == 50.0  # half the bases meet min_depth


def test_quantify_min_depth_boundary_is_inclusive():
    bam = MockBam({"contig1": 100}, default_depth=3)
    ranges = {"contig1": [(10, 20, "geneA")]}
    # depth exactly at the threshold counts as covered
    _, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=3)
    assert coverage["geneA"] == 100.0
    # one above the observed depth counts as uncovered
    _, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=4)
    assert coverage["geneA"] == 0.0


def test_quantify_passes_min_quality_through_to_count_coverage():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneA")]}
    gene_coverage.quantify_gene_coverage(bam, ranges, min_quality=30)
    assert bam.seen_quality == [30]


def test_quantify_sorts_output_by_query():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneB"), (30, 40, "geneA")]}
    depth, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert list(depth) == ["geneA", "geneB"]


def test_quantify_returns_empty_for_no_ranges():
    bam = MockBam({"contig1": 100}, default_depth=5)
    depth, coverage = gene_coverage.quantify_gene_coverage(bam, {}, min_depth=1)
    assert depth == {}
    assert coverage == {}


@pytest.mark.parametrize(
    "contig, segment, message_fragment",
    [
        ("contig1", (10, 10), "start (10) must be < end (10)"),
        ("contig1", (5, 15), "end (15) exceeds contig length (10)"),
        ("missing_contig", (1, 5), "not found in BAM references"),
    ],
)
def test_quantify_edge_guards_raise_clean_value_errors(contig, segment, message_fragment):
    bam = MockBam({"contig1": 10}, default_depth=5)
    ranges = {contig: [(segment[0], segment[1], "geneA")]}
    with pytest.raises(ValueError, match=re.escape(message_fragment)):
        gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)


def test_quantify_raises_on_negative_start():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(-1, 5, "geneA")]}
    with pytest.raises(ValueError, match=re.escape("start (-1) must be >= 0")):
        gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)


# --------------------------------------------------------------------------- #
# integration: GFF parse -> query match -> flatten -> quantify
# --------------------------------------------------------------------------- #

def _write_gff(tmp_path):
    """A two-gene GFF whose CDS carry ``product`` names and one ``exon``."""
    gff = tmp_path / "hier.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\t.\tgene\t1\t40\t.\t+\t.\tID=gene-A\n"
        "contig1\t.\tmRNA\t1\t40\t.\t+\t.\tID=rna-A;Parent=gene-A;product=geneA\n"
        "contig1\t.\tCDS\t11\t20\t.\t+\t0\tID=cds-A;Parent=rna-A;product=geneA\n"
        "contig1\t.\texon\t11\t20\t.\t+\t.\tID=exon-A;Parent=rna-A;product=geneA\n"
        "contig1\t.\tgene\t50\t90\t.\t+\t.\tID=gene-B\n"
        "contig1\t.\tmRNA\t50\t90\t.\t+\t.\tID=rna-B;Parent=gene-B;product=geneB\n"
        "contig1\t.\tCDS\t61\t70\t.\t+\t0\tID=cds-B;Parent=rna-B;product=geneB\n"
    )
    return str(gff)


def test_gff_query_ranges_matches_product_and_quantifies(tmp_path):
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    features = assimilate_gff(_write_gff(tmp_path))
    # match the CDS product; keep only geneA, label by the matched query term
    ranges = gff_query_ranges(
        features, ["geneA"], "RNA", "CDS", ["product"], exact_match=False
    )
    # CDS 1-based [11, 20] -> 0-based half-open [10, 20)
    assert ranges == {"contig1": [(10, 20, "geneA")]}

    bam = MockBam({"contig1": 100}, default_depth=7)
    depth, coverage = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    # geneB was filtered out; only the matched query is reported
    assert depth == {"geneA": 7.0}
    assert coverage == {"geneA": 100.0}


def test_gff_query_ranges_honors_feature_type(tmp_path):
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    features = assimilate_gff(_write_gff(tmp_path))
    # geneA carries both a CDS and an exon over [10, 20); selecting exon still works
    exon_ranges = gff_query_ranges(
        features, ["geneA"], "RNA", "exon", ["product"], exact_match=False
    )
    assert exon_ranges == {"contig1": [(10, 20, "geneA")]}

    # a feature_type absent from the matched unit resolves no coordinates
    empty = gff_query_ranges(
        features, ["geneB"], "RNA", "exon", ["product"], exact_match=False
    )
    assert empty == {}


def test_quantify_with_real_bam_counts_a_single_read(make_bam):
    from theiagene.lib.parsers import import_bam

    bam_path = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    imported = import_bam(bam_path)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage = gene_coverage.quantify_gene_coverage(imported, ranges, min_depth=1)

    # the single 50 bp read fully spans [10, 20) at depth 1
    assert depth["geneA"] == 1.0
    assert coverage["geneA"] == 100.0
