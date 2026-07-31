"""Unit tests for theiagene.lib.query (query-name resolution helpers)."""

import pytest

from theiagene.lib import query


def test_exact_and_substring_check():
    assert query.exact_check({"geneA", "geneB"}, "geneA")
    assert not query.exact_check({"geneA"}, "geneAB")
    assert query.substring_check({"geneA"}, "prefix_geneA_suffix")
    assert not query.substring_check({"geneA"}, "geneB")


def test_extract_queries_from_bed(tmp_path):
    bed = tmp_path / "q.bed"
    bed.write_text("contig1\t0\t10\tgeneA\ncontig1\t20\t30\tgeneB\n")
    assert query.extract_queries_from_bed(str(bed)) == {"geneA", "geneB"}


def test_normalize_name_collapses_reserved_characters():
    assert (
        query.normalize_name("lanosterol 14-alpha demethylase")
        == "lanosterol.14-alpha.demethylase"
    )
    # adjacent reserved chars collapse to a single dot, edges stripped
    assert query.normalize_name("a;;b==c") == "a.b.c"
    assert query.normalize_name("  spaced  ") == "spaced"
    assert query.normalize_name(",leading,and,trailing,") == "leading.and.trailing"


def test_sanitize_info_value_replaces_reserved_characters():
    # VCF INFO fields cannot carry whitespace/';'/'='/','; each becomes '_'
    assert query.sanitize_info_value("a b;c=d,e") == "a_b_c_d_e"


def test_match_query_is_normalization_aware():
    # a dotted query matches a spaced identifier and vice versa
    assert (
        query.match_query(
            ["lanosterol.14-alpha.demethylase"],
            ["lanosterol 14-alpha demethylase"],
            exact_match=True,
        )
        == "lanosterol.14-alpha.demethylase"
    )
    # substring match returns the query term that matched
    assert query.match_query(["FKS"], ["FKS1"], exact_match=False) == "FKS"
    # an exact query does not substring-leak under exact_match
    assert query.match_query(["FKS"], ["FKS1"], exact_match=True) is None
    assert query.match_query(["FKS"], ["ERG11"], exact_match=False) is None


def test_match_query_returns_first_query_in_order():
    # ordering follows the query list, not the identifier list
    assert query.match_query(["b", "a"], ["a", "b"], exact_match=True) == "b"


def test_ordered_query_genes_flattens_and_dedupes():
    assert query.ordered_query_genes(["a,b", "b,c", " d "]) == ["a", "b", "c", "d"]
    assert query.ordered_query_genes(None) == []
    assert query.ordered_query_genes([]) == []
    # blank chunks are dropped
    assert query.ordered_query_genes(["a,,b", " "]) == ["a", "b"]
