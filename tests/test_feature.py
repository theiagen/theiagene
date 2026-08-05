"""Unit tests for theiagene.lib.feature (the shared Feature data model)."""

import pytest

from theiagene.lib.feature import (
    Feature,
    FeatureCol,
    group_features,
    _as_int,
    _as_strand,
    _escape_gff_attr,
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
    with pytest.raises(ValueError, match="start is not an integer"):
        _as_int("abc", "start")


@pytest.mark.parametrize(
    "value, expected",
    [("+", 1), ("-", -1), (1, 1), (-1, -1), (".", None), ("?", None), ("", None), (None, None)],
)
def test_as_strand_coerces_all_accepted_forms(value, expected):
    assert _as_strand(value) == expected


def test_as_strand_rejects_unknown_token():
    with pytest.raises(ValueError, match="strand is not a valid GFF3 strand"):
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


@pytest.mark.parametrize(
    "start, end",
    [(None, 6), (0, None), (None, None), (".", 6), (0, "")],
)
def test_feature_requires_both_coordinates(start, end):
    # a missing/undefined start or end cannot yield a coordinate span
    with pytest.raises(ValueError, match="start and stop coordinates required"):
        Feature(fid="a", start=start, end=end)


def test_feature_accepts_zero_coordinate():
    # start=0 is a valid 0-based coordinate, not a missing value
    feat = Feature(fid="a", start=0, end=6)
    assert (feat.start, feat.end) == (0, 6)


def test_feature_defaults_are_not_shared_between_instances():
    # attributes/descendants must be per-instance, else the hierarchy collapses
    a = Feature(fid="a", start=0, end=6)
    b = Feature(fid="b", start=0, end=6)
    a.attributes["k"] = "v"
    a.descendants.append(b)
    assert b.attributes == {}
    assert b.descendants == []


# --------------------------------------------------------------------------- #
# attribute serialization
# --------------------------------------------------------------------------- #

def test_gff_escape_percent_encodes_reserved_first():
    # '%' is encoded first so later escapes are not themselves re-encoded
    assert _escape_gff_attr("a=b;c,d") == "a%3Db%3Bc%2Cd"
    assert _escape_gff_attr("100%") == "100%25"


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
    # group_features keys off the ID/Parent attributes (as real GFF records carry
    # them), so mirror that here rather than relying on the fid/pid arguments alone
    attributes = {"ID": fid}
    if pid:
        attributes["Parent"] = pid
    return Feature(fid=fid, pid=pid, seqid=seqid, type=type, start=0, end=6,
                   attributes=attributes)


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


def test_group_features_dedupes_duplicate_ids():
    # a repeated id with nobody referencing it as a parent is made unique, not rejected
    grouped = group_features([_feat("dup", "gene"), _feat("dup", "gene")])
    assert {f.fid for f in grouped["chr1"]} == {"dup", "dup_1"}


def test_group_features_rejects_ambiguous_parent_id():
    # 'dup' names two features *and* is referenced as a Parent -> the link is ambiguous
    feats = [_feat("dup", "gene"), _feat("dup", "gene"), _feat("kid", "CDS", pid="dup")]
    with pytest.raises(KeyError, match="belongs to multiple features"):
        group_features(feats)


def test_group_features_falls_back_to_fid_without_id_attribute():
    # features built directly (explicit fid=, no ID attribute) must keep their
    # own identity: without the fallback every one resolves to None, collides on
    # the None key, and is silently renamed None_1, None_2, ...
    a = Feature(seqid="c1", start=0, end=100, type="gene", fid="geneA")
    b = Feature(seqid="c1", start=200, end=300, type="gene", fid="geneB")
    grouped = group_features([a, b])
    assert a.fid == "geneA" and b.fid == "geneB"
    assert [f.fid for f in grouped["c1"]] == ["geneA", "geneB"]


def test_featurecol_by_id_resolves_directly_built_features():
    # regression: FeatureCol([...]) from a plain list[Feature] indexes each by its
    # supplied fid; before the fix by_id raised KeyError for every feature past the
    # first because their fids had been overwritten to None_1, None_2, ...
    a = Feature(seqid="c1", start=0, end=100, type="gene", fid="geneA")
    b = Feature(seqid="c1", start=200, end=300, type="gene", fid="geneB")
    col = FeatureCol([a, b])
    assert col.by_id("geneA") is a
    assert col.by_id("geneB") is b


def test_featurecol_buckets_features_by_canonical_class():
    feats = [_feat("a", "mRNA"), _feat("b", "CDS"), _feat("c", "ncRNA")]
    # unrelated features (no Parent/ID links); skip grouping
    fl = FeatureCol(feats, group=False)
    # any *rna type collapses to the 'rna' class, keyed by class name or a raw type
    assert {f.fid for f in fl["rna"]} == {"a", "c"}
    assert {f.fid for f in fl["mRNA"]} == {"a", "c"}
    assert [f.fid for f in fl["cds"]] == ["b"]
    # a class with no members is empty
    assert fl["exon"] == []
