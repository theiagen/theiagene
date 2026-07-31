"""Shared test setup and fixtures.

Puts the src-layout package on ``sys.path`` (so tests run without an install)
and provides a small pysam-BAM factory used by the coverage tests."""

import os
import sys

import pytest

SRC = os.path.join(os.path.dirname(__file__), "src")
if SRC not in sys.path:
    sys.path.insert(0, SRC)


@pytest.fixture
def make_bam(tmp_path):
    """Return a factory that writes an indexed single-read BAM and its path.

    The read is a perfect (all-match, high-quality) alignment, so every base it
    spans has depth 1 -- enough to exercise real ``pysam.count_coverage``."""
    import pysam

    def _make(name="reads.bam", contig="contig1", contig_len=100,
              read_start=10, read_len=50):
        path = tmp_path / name
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
        return str(path)

    return _make
