"""Unit tests for theiagene.lib.feature (the shared Feature data model)."""

import pytest

from theiagene.lib.feature import (
    Feature,
    group_features,
    extract_features,
    _as_int,
    _as_strand,
    _gff_escape,
    _format_gff_attributes,
)


# --------------------------------------------------------------------------- #
# column coercion helpers
# --------------------------------------------------------------------------- #

def test_as_int_maps_undefined_to_none_and_parses_digits():
    assert _as_int(".", "start") is None
    assert _as_int("", "start") is None
    assert _as_int(None, "start") is None
    assert _as_int("41", "start") == 41
    assert _as_int(41, "start") == 41


def test_as_int_raises_on_non_numeric():
    with pytest.raises(ValueError, match="start must be an integer"):
        _as_int("abc", "start")


@pytest.mark.parametrize(
    "value, expected",
    [("+", 1), ("-", -1), (1, 1), (-1, -1), (".", None), ("?", None), ("", None), (None, None)],
)
def test_as_strand_coerces_all_accepted_forms(value, expected):
    assert _as_strand(value) == expected


def test_as_strand_rejects_unknown_token():
    with pytest.raises(ValueError, match="strand must be one of"):
        _as_strand("x")


# --------------------------------------------------------------------------- #
# Feature construction
# --------------------------------------------------------------------------- #

def test_feature_sorts_coordinates_and_coerces_columns():
    # start/end are stored smallest-first regardless of input order
    feat = Feature(
        fid="cds-A", pid="rna-A", seqid="chr1", source=".", type="CDS",
        start=20, end=9, score=".", strand="-", phase="0",
    )
    assert (feat.start, feat.end) == (9, 20)
    assert feat.strand == -1
    assert feat.score is None      # '.' -> None
    assert feat.phase == 0
    assert feat.fid == "cds-A" and feat.pid == "rna-A"


def test_feature_defaults_are_not_shared_between_instances():
    # attributes/descendants must be per-instance, else the hierarchy collapses
    a = Feature(fid="a")
    b = Feature(fid="b")
    a.attributes["k"] = "v"
    a.descendants.append(b)
    assert b.attributes == {}
    assert b.descendants == []


# --------------------------------------------------------------------------- #
# attribute serialization
# --------------------------------------------------------------------------- #

def test_gff_escape_percent_encodes_reserved_first():
    # '%' is encoded first so later escapes are not themselves re-encoded
    assert _gff_escape("a=b;c,d") == "a%3Db%3Bc%2Cd"
    assert _gff_escape("100%") == "100%25"


def test_format_gff_attributes_roundtrips_and_handles_empty():
    assert _format_gff_attributes({}) == "."
    assert _format_gff_attributes({"ID": "x", "product": "y"}) == "ID=x;product=y"


# --------------------------------------------------------------------------- #
# GFF3 serialization (to_gff)
# --------------------------------------------------------------------------- #

def test_to_gff_line_restores_one_based_inclusive_coordinates():
    # internal 0-based half-open [9, 33) -> GFF 1-based inclusive [10, 33]
    feat = Feature(
        fid="cds-A", pid="rna-A", seqid="chr1", source=".", type="CDS",
        start=9, end=33, score=None, strand="+", phase=0,
        attributes={"product": "alpha"},
    )
    assert feat.to_gff() == (
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=cds-A;Parent=rna-A;product=alpha"
    )


def test_to_gff_emits_descendants_depth_first():
    parent = Feature(
        fid="rna-A", seqid="chr1", source=".", type="mRNA",
        start=0, end=6, strand="+", phase=".",
    )
    child = Feature(
        fid="cds-A", pid="rna-A", seqid="chr1", source=".", type="CDS",
        start=0, end=6, strand="+", phase=0,
    )
    parent.descendants.append(child)
    lines = parent.to_gff().splitlines()
    assert lines[0].startswith("chr1\t.\tmRNA")
    assert lines[1].startswith("chr1\t.\tCDS")
    assert "Parent=rna-A" in lines[1]


# --------------------------------------------------------------------------- #
# grouping and filtering
# --------------------------------------------------------------------------- #

def _feat(fid, type, pid=None, seqid="chr1"):
    return Feature(fid=fid, pid=pid, seqid=seqid, type=type, start=0, end=6)


def test_group_features_wires_parent_and_descendant_links():
    gene = _feat("gene-A", "gene")
    rna = _feat("rna-A", "mRNA", pid="gene-A")
    cds = _feat("cds-A", "CDS", pid="rna-A")
    grouped = group_features([gene, rna, cds])

    # every feature on the contig is returned under its seqid
    assert set(f.fid for f in grouped["chr1"]) == {"gene-A", "rna-A", "cds-A"}
    # parent/descendant links follow the Parent chain
    assert rna.parent is gene and cds.parent is rna
    assert rna.descendants == [cds]
    assert gene.descendants == [rna]


def test_group_features_keys_contigs_separately():
    a = _feat("gene-A", "gene", seqid="chr1")
    b = _feat("gene-B", "gene", seqid="chr2")
    grouped = group_features([a, b])
    assert set(grouped) == {"chr1", "chr2"}
    assert [f.fid for f in grouped["chr2"]] == ["gene-B"]


def test_group_features_rejects_duplicate_ids():
    with pytest.raises(KeyError, match="dup"):
        group_features([_feat("dup", "gene"), _feat("dup", "gene")])


def test_extract_features_substring_vs_exact():
    feats = [_feat("a", "mRNA"), _feat("b", "CDS"), _feat("c", "ncRNA")]
    # substring: "RNA" matches both mRNA and ncRNA
    assert {f.fid for f in extract_features(feats, "RNA")} == {"a", "c"}
    # exact: only a literal type match
    assert [f.fid for f in extract_features(feats, "mRNA", exact_match=True)] == ["a"]
    assert extract_features(feats, "exon") == []
