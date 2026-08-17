"""Unit tests for theiagene.gene_coverage (coverage quantification + helpers).

Coverage is quantified over a ``{<CONTIG>: [(START, END, LABEL), ...]}`` map (the
shape :func:`theiagene.lib.query.gff_query_ranges`/``bed_query_ranges`` produce);
every range sharing a ``LABEL`` is pooled onto one output row."""

import argparse
import json
import re

import pytest

from theiagene import gene_coverage


# --------------------------------------------------------------------------- #
# test doubles
# --------------------------------------------------------------------------- #

class MockRead:
    """A minimal stand-in for pysam.AlignedSegment (what ``fetch`` yields)."""

    def __init__(self, query_name, contig, start, length=10, mapping_quality=60,
                 flag=0, is_unmapped=False, is_secondary=False,
                 is_supplementary=False, is_qcfail=False, is_duplicate=False):
        self.query_name = query_name
        self.reference_name = contig
        self.reference_start = start
        self.reference_end = start + length
        self.mapping_quality = mapping_quality
        self.flag = flag
        self.is_unmapped = is_unmapped
        self.is_secondary = is_secondary
        self.is_supplementary = is_supplementary
        self.is_qcfail = is_qcfail
        self.is_duplicate = is_duplicate


class MockBam:
    """A minimal stand-in for pysam.AlignmentFile.

    ``depth_fn(contig, pos) -> int`` supplies per-base depth (default: a flat
    ``default_depth``); all of it is reported in the first (A) channel so the
    four-channel sum used by quantify_gene_coverage equals that depth.
    ``reads`` are MockReads returned by ``fetch`` when they overlap the queried
    region -- depth and reads-mapped are independent here, so a test can set
    either without the other."""

    def __init__(self, lengths, depth_fn=None, default_depth=5, reads=None):
        self._lengths = dict(lengths)
        self.references = tuple(lengths)
        self._depth_fn = depth_fn
        self._default = default_depth
        self._reads = list(reads or [])
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

    def fetch(self, contig, start, end):
        for read in self._reads:
            if (read.reference_name == contig
                    and read.reference_start < end and read.reference_end > start):
                yield read


class MockMultiChannelBam:
    """A MockBam whose per-base depth is spread across all four count_coverage
    channels, so a test can prove quantify sums every channel (not just A)."""

    def __init__(self, lengths, per_channel=(1, 2, 3, 4)):
        self._lengths = dict(lengths)
        self.references = tuple(lengths)
        self._per_channel = per_channel

    def get_reference_length(self, contig):
        return self._lengths[contig]

    def count_coverage(self, contig, start, end, quality_threshold=0):
        width = end - start
        return tuple([value] * width for value in self._per_channel)

    def fetch(self, contig, start, end):
        return ()


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
    with pytest.raises(FileNotFoundError, match="query_genes"):
        gene_coverage.input_error_handling(_Args(reference_gff="ref.gff"))


# --------------------------------------------------------------------------- #
# make_tsv
# --------------------------------------------------------------------------- #

def test_make_tsv_renders_header_and_rows():
    tsv = gene_coverage.make_tsv(
        {"geneA": 5.0, "geneB": 2.0},
        {"geneA": 100.0, "geneB": 50.0},
        {"geneA": 12, "geneB": 3},
        {"geneA": True, "geneB": False},
        ambiguous_contig=False,
    )
    lines = tsv.splitlines()
    assert lines[0] == (
        "#query\taverage_depth\tpercent_coverage\treads_mapped\treads_mapped_pass"
    )
    assert lines[1] == "geneA\t5.0\t100.0\t12\tTrue"
    assert lines[2] == "geneB\t2.0\t50.0\t3\tFalse"


def test_make_tsv_warns_in_header_when_contig_is_ambiguous():
    tsv = gene_coverage.make_tsv(
        {"geneA": 1.0}, {"geneA": 100.0}, {"geneA": 1}, {"geneA": True},
        ambiguous_contig=True,
    )
    assert tsv.splitlines()[0].startswith("#query (WARNING")


def test_make_tsv_with_no_rows_is_header_only():
    tsv = gene_coverage.make_tsv({}, {}, {}, {}, ambiguous_contig=False)
    assert tsv == (
        "#query\taverage_depth\tpercent_coverage\treads_mapped\treads_mapped_pass\n"
    )


# --------------------------------------------------------------------------- #
# flag_reads_mapped
# --------------------------------------------------------------------------- #

def test_flag_reads_mapped_is_inclusive_at_the_threshold():
    flags = gene_coverage.flag_reads_mapped(
        {"below": 9, "at": 10, "above": 11}, min_reads_mapped=10
    )
    assert flags == {"below": False, "at": True, "above": True}


def test_flag_reads_mapped_defaults_to_one_read():
    # the default separates a query with any evidence from one with none
    assert gene_coverage.flag_reads_mapped({"none": 0, "some": 1}) == {
        "none": False,
        "some": True,
    }


def test_flag_reads_mapped_keeps_every_query():
    # the flag is a call layered onto the measurements, never a filter
    flags = gene_coverage.flag_reads_mapped({"geneA": 100, "geneB": 0}, 5)
    assert set(flags) == {"geneA", "geneB"}


# --------------------------------------------------------------------------- #
# union_intervals
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize(
    "intervals, expected",
    [
        ([(0, 50), (0, 100)], [(0, 100)]),            # overlap -> single span
        ([(0, 50), (50, 100)], [(0, 100)]),           # contiguous -> merged
        ([(0, 50), (60, 100)], [(0, 50), (60, 100)]),  # disjoint -> kept apart
        ([(60, 100), (0, 50)], [(0, 50), (60, 100)]),  # unsorted input -> ordered
        ([(0, 100), (20, 40)], [(0, 100)]),           # nested -> outer span
    ],
)
def test_union_intervals_merges_overlapping_and_contiguous(intervals, expected):
    assert gene_coverage.union_intervals(intervals) == expected


# --------------------------------------------------------------------------- #
# quantify_gene_coverage
# --------------------------------------------------------------------------- #

def test_quantify_known_depth_and_full_breadth():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    assert depth["geneA"] == 5.0
    assert coverage["geneA"] == 100.0


def test_quantify_pools_segments_sharing_a_label():
    bam = MockBam({"contig1": 100}, default_depth=4)
    # two CDS segments of one query gene share a label and are pooled
    ranges = {"contig1": [(0, 6, "geneA"), (10, 16, "geneA")]}

    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    # 12 evaluated bases, all at depth 4
    assert depth["geneA"] == 4.0
    assert coverage["geneA"] == 100.0


def test_quantify_unions_overlapping_ranges_sharing_a_label():
    # regression: two overlapping isoform CDS under one label must be counted
    # over their per-base union, not concatenated (which double-weights the
    # shared bases). Depth 1 over [0, 50), uncovered over [50, 100).
    def depth_fn(contig, pos):
        return 1 if pos < 50 else 0

    bam = MockBam({"c1": 200}, depth_fn=depth_fn)
    ranges = {"c1": [(0, 50, "GENE"), (0, 100, "GENE")]}

    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    # union is [0, 100): 100 bases, 50 covered -> 50% breadth, 0.5x depth.
    # a 150-base concatenation would wrongly report 66.67% / 0.6667x.
    assert coverage["GENE"] == 50.0
    assert depth["GENE"] == 0.5


def test_quantify_partial_breadth_from_min_depth_threshold():
    # first half of the region at depth 10, second half uncovered
    def depth_fn(contig, pos):
        return 10 if pos < 15 else 0

    bam = MockBam({"contig1": 100}, depth_fn=depth_fn)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)

    assert depth["geneA"] == 5.0     # mean of five 10s and five 0s
    assert coverage["geneA"] == 50.0  # half the bases meet min_depth


def test_quantify_min_depth_boundary_is_inclusive():
    bam = MockBam({"contig1": 100}, default_depth=3)
    ranges = {"contig1": [(10, 20, "geneA")]}
    # depth exactly at the threshold counts as covered
    _, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=3)
    assert coverage["geneA"] == 100.0
    # one above the observed depth counts as uncovered
    _, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=4)
    assert coverage["geneA"] == 0.0


def test_quantify_passes_min_quality_through_to_count_coverage():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneA")]}
    gene_coverage.quantify_gene_coverage(bam, ranges, min_quality=30)
    assert bam.seen_quality == [30]


def test_quantify_sorts_output_by_query():
    bam = MockBam({"contig1": 100}, default_depth=5)
    ranges = {"contig1": [(10, 20, "geneB"), (30, 40, "geneA")]}
    depth, _, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert list(depth) == ["geneA", "geneB"]


def test_quantify_returns_empty_for_no_ranges():
    bam = MockBam({"contig1": 100}, default_depth=5)
    depth, coverage, reads = gene_coverage.quantify_gene_coverage(bam, {}, min_depth=1)
    assert depth == {}
    assert coverage == {}
    assert reads == {}


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


def test_quantify_seeds_expected_labels_absent_from_ranges():
    bam = MockBam({"contig1": 100}, default_depth=5,
                  reads=[MockRead("r1", "contig1", 10)])
    ranges = {"contig1": [(10, 20, "geneA")]}
    depth, coverage, reads = gene_coverage.quantify_gene_coverage(
        bam, ranges, min_depth=1, expected_labels=["geneA", "geneMissing"]
    )
    # geneMissing resolved no coordinates but is still reported as absent (0)
    assert depth == {"geneA": 5.0, "geneMissing": 0}
    assert coverage == {"geneA": 100.0, "geneMissing": 0}
    assert reads == {"geneA": 1, "geneMissing": 0}


def test_quantify_measured_label_overrides_the_zero_seed():
    bam = MockBam({"contig1": 100}, default_depth=8)
    ranges = {"contig1": [(10, 20, "geneA")]}
    depth, coverage, _ = gene_coverage.quantify_gene_coverage(
        bam, ranges, min_depth=1, expected_labels=["geneA"]
    )
    # the seeded 0 is replaced by the measured value, not summed with it
    assert depth["geneA"] == 8.0
    assert coverage["geneA"] == 100.0


def test_quantify_all_expected_labels_zero_when_no_ranges():
    bam = MockBam({"contig1": 100}, default_depth=5)
    depth, coverage, reads = gene_coverage.quantify_gene_coverage(
        bam, {}, min_depth=1, expected_labels=["geneB", "geneA"]
    )
    # nothing measured; every requested query is reported as absent and sorted
    assert depth == {"geneA": 0, "geneB": 0}
    assert coverage == {"geneA": 0, "geneB": 0}
    assert reads == {"geneA": 0, "geneB": 0}
    assert list(depth) == ["geneA", "geneB"]


# --------------------------------------------------------------------------- #
# reads mapped
# --------------------------------------------------------------------------- #

def test_quantify_counts_reads_overlapping_a_range():
    bam = MockBam(
        {"contig1": 100},
        reads=[
            MockRead("in_range", "contig1", 12),
            MockRead("elsewhere", "contig1", 60),  # outside [10, 20)
        ],
    )
    ranges = {"contig1": [(10, 20, "geneA")]}
    _, _, reads = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert reads == {"geneA": 1}


def test_quantify_counts_a_read_spanning_pooled_ranges_once():
    # one 40 bp read overlaps both CDS segments of geneA; pooling the label's
    # ranges must not count it twice
    bam = MockBam({"contig1": 100}, reads=[MockRead("spanner", "contig1", 5, length=40)])
    ranges = {"contig1": [(10, 20, "geneA"), (30, 40, "geneA")]}
    _, _, reads = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert reads == {"geneA": 1}


def test_quantify_counts_paired_mates_separately():
    # mates share a query name, so the flag is what keeps them distinct
    bam = MockBam(
        {"contig1": 100},
        reads=[
            MockRead("pair", "contig1", 10, flag=99),
            MockRead("pair", "contig1", 14, flag=147),
        ],
    )
    ranges = {"contig1": [(10, 20, "geneA")]}
    _, _, reads = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert reads == {"geneA": 2}


def test_quantify_min_mapping_quality_filters_reads():
    bam = MockBam(
        {"contig1": 100},
        reads=[
            MockRead("high", "contig1", 12, mapping_quality=60),
            MockRead("low", "contig1", 12, mapping_quality=5),
        ],
    )
    ranges = {"contig1": [(10, 20, "geneA")]}

    _, _, unfiltered = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert unfiltered == {"geneA": 2}

    _, _, filtered = gene_coverage.quantify_gene_coverage(
        bam, ranges, min_depth=1, min_mapping_quality=20
    )
    assert filtered == {"geneA": 1}


def test_quantify_min_mapping_quality_boundary_is_inclusive():
    bam = MockBam({"contig1": 100},
                  reads=[MockRead("edge", "contig1", 12, mapping_quality=30)])
    ranges = {"contig1": [(10, 20, "geneA")]}
    # mapping quality exactly at the threshold still counts
    _, _, at = gene_coverage.quantify_gene_coverage(
        bam, ranges, min_depth=1, min_mapping_quality=30
    )
    assert at == {"geneA": 1}
    # one above the observed mapping quality does not
    _, _, above = gene_coverage.quantify_gene_coverage(
        bam, ranges, min_depth=1, min_mapping_quality=31
    )
    assert above == {"geneA": 0}


@pytest.mark.parametrize(
    "flag_kwarg",
    ["is_unmapped", "is_secondary", "is_supplementary", "is_qcfail", "is_duplicate"],
)
def test_quantify_excludes_reads_count_coverage_ignores(flag_kwarg):
    # the reads total must describe the same alignments that feed depth/breadth
    bam = MockBam(
        {"contig1": 100},
        reads=[
            MockRead("counted", "contig1", 12),
            MockRead("excluded", "contig1", 12, **{flag_kwarg: True}),
        ],
    )
    ranges = {"contig1": [(10, 20, "geneA")]}
    _, _, reads = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert reads == {"geneA": 1}


def test_quantify_pools_reads_for_one_label_across_contigs():
    bam = MockBam(
        {"contig1": 100, "contig2": 100},
        reads=[
            MockRead("r1", "contig1", 0),
            MockRead("r2", "contig2", 0),
        ],
    )
    ranges = {"contig1": [(0, 10, "geneA")], "contig2": [(0, 10, "geneA")]}
    _, _, reads = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert reads == {"geneA": 2}


def test_quantify_across_multiple_contigs_keeps_labels_separate():
    bam = MockBam({"contig1": 100, "contig2": 100}, default_depth=6)
    ranges = {
        "contig1": [(10, 20, "geneA")],
        "contig2": [(30, 40, "geneB")],
    }
    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    assert depth == {"geneA": 6.0, "geneB": 6.0}
    assert coverage == {"geneA": 100.0, "geneB": 100.0}


def test_quantify_pools_one_label_across_contigs():
    # depth 10 over contig1's bases, depth 0 over contig2's bases
    def depth_fn(contig, pos):
        return 10 if contig == "contig1" else 0

    bam = MockBam({"contig1": 100, "contig2": 100}, depth_fn=depth_fn)
    ranges = {
        "contig1": [(0, 10, "geneA")],
        "contig2": [(0, 10, "geneA")],
    }
    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    # both contigs' bases pool onto one geneA row: mean 5, half covered
    assert depth["geneA"] == 5.0
    assert coverage["geneA"] == 50.0


def test_quantify_sums_all_four_coverage_channels():
    bam = MockMultiChannelBam({"contig1": 100}, per_channel=(1, 2, 3, 4))
    ranges = {"contig1": [(10, 20, "geneA")]}
    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    # each base's depth is 1+2+3+4 = 10 summed across the A/C/G/T channels
    assert depth["geneA"] == 10.0
    assert coverage["geneA"] == 100.0


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
    # match the CDS product; keep only geneA, label by its resolved gene name
    ranges = gff_query_ranges(
        features, ["geneA"], "CDS", ["product"], exact_match=False
    )
    # CDS 1-based [11, 20] -> 0-based half-open [10, 20)
    assert ranges == {"contig1": [(10, 20, "geneA")]}

    bam = MockBam({"contig1": 100}, default_depth=7)
    depth, coverage, _ = gene_coverage.quantify_gene_coverage(bam, ranges, min_depth=1)
    # geneB was filtered out; only the matched query is reported
    assert depth == {"geneA": 7.0}
    assert coverage == {"geneA": 100.0}


def _write_paralog_gff(tmp_path):
    """Two paralogs (ERG1, ERG11) whose gene name lives on the parent ``gene``
    record and CDS ``product`` -- not on the mRNA. ``ERG1`` is a substring of
    ``ERG11``, so a substring query for ``ERG1`` matches both units."""
    gff = tmp_path / "paralogs.gff"
    gff.write_text(
        "c1\t.\tgene\t1\t300\t.\t+\t.\tID=g1;gene=ERG1\n"
        "c1\t.\tmRNA\t1\t300\t.\t+\t.\tID=r1;Parent=g1\n"
        "c1\t.\tCDS\t1\t300\t.\t+\t0\tID=c1a;Parent=r1;product=ERG1\n"
        "c1\t.\tgene\t1000\t1300\t.\t+\t.\tID=g2;gene=ERG11\n"
        "c1\t.\tmRNA\t1000\t1300\t.\t+\t.\tID=r2;Parent=g2\n"
        "c1\t.\tCDS\t1000\t1300\t.\t+\t0\tID=c2a;Parent=r2;product=ERG11\n"
    )
    return str(gff)


def test_gff_query_ranges_labels_paralogs_by_resolved_gene_not_query_term(tmp_path):
    # regression: a substring query matching two paralogs must yield two rows
    # labelled by each unit's resolved gene name -- not one row labelled by the
    # shared query term (which would collapse distinct genes together)
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    features = assimilate_gff(_write_paralog_gff(tmp_path))
    ranges = gff_query_ranges(
        features, ["ERG1"], "CDS", ["product"], exact_match=False
    )
    # both mRNAs match the "ERG1" substring, but resolve to distinct labels
    assert ranges == {"c1": [(0, 300, "ERG1"), (999, 1300, "ERG11")]}


def test_gff_query_ranges_falls_back_to_gene_without_rna(tmp_path):
    # regression: Bakta / RefSeq bacterial GFF3 carry no mRNA record -- CDS is a
    # direct child of gene -- so grouping must fall back to the gene level rather
    # than reporting the covered gene as zero coordinates
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    gff = tmp_path / "bakta.gff"
    gff.write_text(
        "c1\t.\tgene\t1\t900\t.\t+\t.\tID=g1;gene=ERG11\n"
        "c1\t.\tCDS\t1\t900\t.\t+\t0\tID=c1cds;Parent=g1;product=ERG11\n"
    )
    ranges = gff_query_ranges(
        assimilate_gff(str(gff)), ["ERG11"], "CDS", ["product"], exact_match=False
    )
    assert ranges == {"c1": [(0, 900, "ERG11")]}


def test_gff_query_ranges_falls_back_to_feature_type_without_gene(tmp_path):
    # an annotation carrying only bare CDS records (no gene parent) resolves by
    # grouping directly on the feature_type
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    gff = tmp_path / "cds_only.gff"
    gff.write_text("c1\t.\tCDS\t1\t900\t.\t+\t0\tID=c1cds;product=ERG11\n")
    ranges = gff_query_ranges(
        assimilate_gff(str(gff)), ["ERG11"], "CDS", ["product"], exact_match=False
    )
    assert ranges == {"c1": [(0, 900, "ERG11")]}


def test_gff_query_ranges_honors_feature_type(tmp_path):
    from theiagene.lib.parsers import assimilate_gff
    from theiagene.lib.query import gff_query_ranges

    features = assimilate_gff(_write_gff(tmp_path))
    # geneA carries both a CDS and an exon over [10, 20); selecting exon still works
    exon_ranges = gff_query_ranges(
        features, ["geneA"], "exon", ["product"], exact_match=False
    )
    assert exon_ranges == {"contig1": [(10, 20, "geneA")]}

    # a feature_type absent from the matched unit resolves no coordinates
    empty = gff_query_ranges(
        features, ["geneB"], "exon", ["product"], exact_match=False
    )
    assert empty == {}


def test_quantify_with_real_bam_counts_a_single_read(make_bam):
    from theiagene.lib.parsers import import_bam

    bam_path = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    imported = import_bam(bam_path)
    ranges = {"contig1": [(10, 20, "geneA")]}

    depth, coverage, reads = gene_coverage.quantify_gene_coverage(
        imported, ranges, min_depth=1
    )

    # the single 50 bp read fully spans [10, 20) at depth 1
    assert depth["geneA"] == 1.0
    assert coverage["geneA"] == 100.0
    assert reads["geneA"] == 1


def test_quantify_with_real_bam_honors_min_mapping_quality(make_bam):
    from theiagene.lib.parsers import import_bam

    # the fixture's read is written with mapping quality 60
    bam_path = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    ranges = {"contig1": [(10, 20, "geneA")]}

    _, _, kept = gene_coverage.quantify_gene_coverage(
        import_bam(bam_path), ranges, min_depth=1, min_mapping_quality=60
    )
    assert kept["geneA"] == 1

    _, _, dropped = gene_coverage.quantify_gene_coverage(
        import_bam(bam_path), ranges, min_depth=1, min_mapping_quality=61
    )
    assert dropped["geneA"] == 0


# --------------------------------------------------------------------------- #
# run_cli / main: end-to-end orchestration and output files
# --------------------------------------------------------------------------- #

def _cli_args(bam, **overrides):
    """A fully-defaulted gene_coverage Namespace with per-test overrides."""
    parser = argparse.ArgumentParser()
    gene_coverage.add_arguments(parser)
    args = parser.parse_args(["--bam", bam])
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


def _read_outputs():
    """Load the five files run_cli writes into the current directory."""
    with open("DEPTH_DICT.json") as fh:
        depth = json.load(fh)
    with open("COVERAGE_DICT.json") as fh:
        coverage = json.load(fh)
    with open("READS_DICT.json") as fh:
        reads = json.load(fh)
    with open("COVERAGE_STATS.tsv") as fh:
        tsv = fh.read()
    return depth, coverage, reads, tsv


def _read_reads_pass():
    """Load the per-query --min_reads_mapped calls run_cli writes."""
    with open("READS_PASS_DICT.json") as fh:
        return json.load(fh)


def _write_bed(tmp_path, rows, name="regions.bed"):
    bed = tmp_path / name
    bed.write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in rows))
    return str(bed)


def _make_multicontig_bam(path, contigs):
    """Write and index a (read-free) BAM declaring several reference contigs."""
    import pysam

    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"LN": length, "SN": name} for name, length in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass
    pysam.index(str(path))
    return str(path)


def test_run_cli_gff_with_query_genes_writes_all_outputs(tmp_path, make_bam, monkeypatch):
    gff = _write_gff(tmp_path)
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(bam, reference_gff=gff, query_genes=["geneA"])

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, coverage, reads, tsv = _read_outputs()
    assert depth == {"geneA": 1.0}
    assert coverage == {"geneA": 100.0}
    assert reads == {"geneA": 1}
    assert _read_reads_pass() == {"geneA": True}
    assert tsv.splitlines()[0] == (
        "#query\taverage_depth\tpercent_coverage\treads_mapped\treads_mapped_pass"
    )
    assert "geneA\t1.0\t100.0\t1\tTrue" in tsv


def test_run_cli_bed_as_coordinate_source(tmp_path, make_bam, monkeypatch):
    # BED is itself the coordinate source: no GFF and no --query_genes
    bed = _write_bed(tmp_path, [("contig1", 10, 20, "geneA")])
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(bam, bedfile=bed)

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, coverage, _, _ = _read_outputs()
    assert depth == {"geneA": 1.0}
    assert coverage == {"geneA": 100.0}


def test_run_cli_gff_source_with_bed_supplied_names(tmp_path, make_bam, monkeypatch):
    # coordinates come from the GFF; the BED only supplies the query names
    gff = _write_gff(tmp_path)
    bed = _write_bed(tmp_path, [("ignored", 0, 1, "geneA")])
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(bam, reference_gff=gff, bedfile=bed)

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, _, _, _ = _read_outputs()
    # geneA resolved from the GFF CDS; the BED coordinates are never consulted
    assert depth == {"geneA": 1.0}


def test_run_cli_ambiguous_contig_collapses_and_warns(tmp_path, make_bam, monkeypatch):
    # BED coordinates sit on a contig whose name differs from the BAM's single
    # contig; ambiguous_contig re-files them onto the BAM contig
    bed = _write_bed(tmp_path, [("ref_chrom", 10, 20, "geneA")])
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(bam, bedfile=bed, ambiguous_contig=True)

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, coverage, reads, tsv = _read_outputs()
    assert depth == {"geneA": 1.0}
    assert coverage == {"geneA": 100.0}
    assert reads == {"geneA": 1}
    # the ambiguous-contig warning is surfaced in the TSV header
    assert tsv.splitlines()[0].startswith("#query (WARNING")


def test_run_cli_ambiguous_contig_rejects_multi_contig_bam(tmp_path, monkeypatch):
    bed = _write_bed(tmp_path, [("ref_chrom", 10, 20, "geneA")])
    bam = _make_multicontig_bam(
        tmp_path / "multi.bam", [("contig1", 100), ("contig2", 100)]
    )
    args = _cli_args(bam, bedfile=bed, ambiguous_contig=True)

    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="multiple contigs"):
        gene_coverage.run_cli(args)


def test_run_cli_reports_unresolved_query_as_zero_and_warns(
    tmp_path, make_bam, monkeypatch, caplog
):
    gff = _write_gff(tmp_path)
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    # a query term that matches no GFF feature resolves no coordinates
    args = _cli_args(bam, reference_gff=gff, query_genes=["ghost_gene"])

    monkeypatch.chdir(tmp_path)
    with caplog.at_level("WARNING"):
        assert gene_coverage.run_cli(args) == 0

    depth, coverage, reads, _ = _read_outputs()
    # still reported (seeded to 0), not silently dropped
    assert depth == {"ghost_gene": 0}
    assert coverage == {"ghost_gene": 0}
    assert reads == {"ghost_gene": 0}
    assert _read_reads_pass() == {"ghost_gene": False}
    assert any("No query-gene coordinates" in r.message for r in caplog.records)


def test_run_cli_min_mapping_quality_reaches_the_reads_count(
    tmp_path, make_bam, monkeypatch
):
    gff = _write_gff(tmp_path)
    # the fixture read has mapping quality 60, below this threshold
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(
        bam, reference_gff=gff, query_genes=["geneA"], min_mapping_quality=61
    )

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, _, reads, _ = _read_outputs()
    # the read is still counted toward depth (min_quality gates bases, not reads)
    assert depth == {"geneA": 1.0}
    assert reads == {"geneA": 0}


def test_run_cli_min_reads_mapped_flags_without_filtering(
    tmp_path, make_bam, monkeypatch
):
    gff = _write_gff(tmp_path)
    # the BAM carries a single read, below this threshold
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)
    args = _cli_args(
        bam, reference_gff=gff, query_genes=["geneA"], min_reads_mapped=10
    )

    monkeypatch.chdir(tmp_path)
    assert gene_coverage.run_cli(args) == 0

    depth, coverage, reads, tsv = _read_outputs()
    # measurements are untouched; only the added flag reports the failure
    assert depth == {"geneA": 1.0}
    assert coverage == {"geneA": 100.0}
    assert reads == {"geneA": 1}
    assert _read_reads_pass() == {"geneA": False}
    assert "geneA\t1.0\t100.0\t1\tFalse" in tsv


def test_main_runs_end_to_end_via_argv(tmp_path, make_bam, monkeypatch):
    gff = _write_gff(tmp_path)
    bam = make_bam(contig="contig1", contig_len=100, read_start=10, read_len=50)

    monkeypatch.chdir(tmp_path)
    rc = gene_coverage.main(
        ["--bam", bam, "--reference_gff", gff, "--query_genes", "geneA"]
    )
    assert rc == 0

    depth, _, _, _ = _read_outputs()
    assert depth == {"geneA": 1.0}
