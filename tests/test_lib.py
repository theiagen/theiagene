"""Unit tests for the shared theiagene.lib libraries."""

import json

import pytest

from theiagene.lib import query, sequence, io_utils, gff


# --------------------------------------------------------------------------- #
# theiagene.lib.query
# --------------------------------------------------------------------------- #

def test_exact_and_substring_check():
    assert query.exact_check({"geneA", "geneB"}, "geneA")
    assert not query.exact_check({"geneA"}, "geneAB")
    assert query.substring_check({"geneA"}, "prefix_geneA_suffix")
    assert not query.substring_check({"geneA"}, "geneB")


def test_normalize_name_collapses_reserved_characters():
    assert (
        query.normalize_name("lanosterol 14-alpha demethylase")
        == "lanosterol.14-alpha.demethylase"
    )
    assert query.normalize_name("a;;b==c") == "a.b.c"
    assert query.normalize_name("  spaced  ") == "spaced"


def test_sanitize_info_value_replaces_reserved_characters():
    assert query.sanitize_info_value("a b;c=d,e") == "a_b_c_d_e"


def test_match_query_is_normalization_aware():
    # dotted query matches a spaced identifier and vice versa
    assert query.match_query(
        ["lanosterol.14-alpha.demethylase"],
        ["lanosterol 14-alpha demethylase"],
        exact_match=True,
    ) == "lanosterol.14-alpha.demethylase"
    assert query.match_query(["FKS"], ["FKS1"], exact_match=False) == "FKS"
    assert query.match_query(["FKS"], ["ERG11"], exact_match=False) is None


def test_ordered_query_genes_flattens_and_dedupes():
    assert query.ordered_query_genes(["a,b", "b,c", " d "]) == ["a", "b", "c", "d"]
    assert query.ordered_query_genes(None) == []


def test_extract_queries_from_bed(tmp_path):
    bed = tmp_path / "q.bed"
    bed.write_text("contig1\t0\t10\tgeneA\ncontig1\t20\t30\tgeneB\n")
    assert query.extract_queries_from_bed(str(bed)) == {"geneA", "geneB"}


# --------------------------------------------------------------------------- #
# theiagene.lib.sequence
# --------------------------------------------------------------------------- #

def test_complement_is_not_reversed():
    assert sequence.complement("ACGT") == "TGCA"
    assert sequence.complement("acgt") == "tgca"


def test_translate_truncates_to_whole_codons():
    assert sequence.translate("ATGTAT", 1) == "MY"
    assert sequence.translate("ATGTA", 1) == "M"  # trailing partial codon dropped
    assert sequence.translate("", 1) == ""


def test_aa3_maps_to_three_letter_with_hgvs_extensions():
    assert sequence.aa3("M") == "Met"
    assert sequence.aa3("*") == "Ter"
    assert sequence.aa3("?") == "Xaa"  # unknown symbol


@pytest.mark.parametrize(
    "allele, expected",
    [
        ("ACGT", True),
        ("acgtn", True),
        ("*", False),        # spanning deletion
        ("<DEL>", False),    # symbolic
        ("", False),         # empty
        (None, False),       # missing
    ],
)
def test_is_nucleotide_allele(allele, expected):
    assert sequence.is_nucleotide_allele(allele) is expected


# --------------------------------------------------------------------------- #
# theiagene.lib.io_utils
# --------------------------------------------------------------------------- #

def test_write_json_writes_data(tmp_path):
    out = tmp_path / "d.json"
    io_utils.write_json(str(out), {"geneA": 5.0})
    assert json.loads(out.read_text()) == {"geneA": 5.0}


def test_write_json_spoofs_cromwell_on_empty(tmp_path):
    out = tmp_path / "empty.json"
    io_utils.write_json(str(out), {})
    assert out.read_text() == '{"": 0}'


# --------------------------------------------------------------------------- #
# theiagene.lib.gff
# --------------------------------------------------------------------------- #

def test_parse_gff_attributes_splits_and_percent_decodes():
    parsed = gff.parse_gff_attributes("ID=cds-1;product=lanosterol%2014-alpha;gene=ERG11;")
    assert parsed == {
        "ID": "cds-1",
        "product": "lanosterol 14-alpha",  # %20 decoded to a space
        "gene": "ERG11",
    }


def test_parse_gff_attributes_ignores_malformed_fields():
    assert gff.parse_gff_attributes("") == {}
    assert gff.parse_gff_attributes("novalue;ID=x") == {"ID": "x"}


def test_iter_gff_features_converts_coordinates_and_strand(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a;product=alpha\n"
        "chr2\t.\tCDS\t6\t17\t.\t-\t0\tID=b;product=beta\n"
    )
    feats = list(gff.iter_gff_features(str(path), "CDS"))
    assert len(feats) == 2
    # GFF 1-based-inclusive [10, 33] -> 0-based half-open [9, 33)
    assert (feats[0]["seqid"], feats[0]["start"], feats[0]["end"], feats[0]["strand"]) == (
        "chr1", 9, 33, 1
    )
    assert feats[1]["strand"] == -1
    assert feats[0]["attributes"]["product"] == "alpha"


def test_iter_gff_features_skips_comments_blanks_fasta_and_other_types(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text(
        "##gff-version 3\n"
        "\n"
        "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=g\n"        # wrong feature type
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a;product=alpha\n"
        "##FASTA\n"
        "chr1\t.\tCDS\t1\t3\t.\t+\t0\tID=ignored\n"     # after ##FASTA: ignored
    )
    feats = list(gff.iter_gff_features(str(path), "CDS"))
    assert len(feats) == 1
    assert feats[0]["attributes"]["ID"] == "a"


def test_iter_gff_features_marks_unresolved_strand_as_none(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text("chr1\t.\tCDS\t10\t33\t.\t.\t0\tID=a\n")
    feats = list(gff.iter_gff_features(str(path), "CDS"))
    assert feats[0]["strand"] is None
