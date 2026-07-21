"""Unit tests for the shared theiagene.lib libraries."""

import json

import pytest

from theiagene.lib import query, sequence, io_utils, gene_model, parsers


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
# theiagene.lib.parsers -- low-level GFF3 readers
# --------------------------------------------------------------------------- #

def test_parse_gff_attributes_splits_and_percent_decodes():
    parsed = parsers.parse_gff_attributes("ID=cds-1;product=lanosterol%2014-alpha;gene=ERG11;")
    assert parsed == {
        "ID": "cds-1",
        "product": "lanosterol 14-alpha",  # %20 decoded to a space
        "gene": "ERG11",
    }


def test_parse_gff_attributes_empty_column_is_empty_dict():
    # an empty (or bare-';') attribute column is valid GFF3, not malformed
    assert parsers.parse_gff_attributes("") == {}
    assert parsers.parse_gff_attributes(";") == {}


def test_parse_gff_attributes_raises_on_malformed_field():
    # a field lacking the value delimiter is malformed and must raise
    with pytest.raises(ValueError, match="unexpected attributes field"):
        parsers.parse_gff_attributes("novalue;ID=x")


def test_parse_gff_attributes_keeps_value_containing_delimiter():
    # only the first '=' splits key from value, so an '=' in the value survives
    assert parsers.parse_gff_attributes("note=a=b") == {"note": "a=b"}


def test_iter_gff_features_converts_coordinates_and_strand(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a;product=alpha\n"
        "chr2\t.\tCDS\t6\t17\t.\t-\t0\tID=b;product=beta\n"
    )
    feats = list(parsers.iter_gff_features(str(path), "CDS"))
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
    feats = list(parsers.iter_gff_features(str(path), "CDS"))
    assert len(feats) == 1
    assert feats[0]["attributes"]["ID"] == "a"


def test_iter_gff_features_marks_unresolved_strand_as_none(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text("chr1\t.\tCDS\t10\t33\t.\t.\t0\tID=a\n")
    feats = list(parsers.iter_gff_features(str(path), "CDS"))
    assert feats[0]["strand"] is None


# --------------------------------------------------------------------------- #
# theiagene.lib.gene_model
# --------------------------------------------------------------------------- #

def test_gene_span_and_positions_are_derived_from_parts():
    plus = gene_model.Gene("g", "chr1", strand=1, parts=[(9, 12), (3, 6)])
    assert plus.genomic_start == 3
    assert plus.genomic_end == 12
    # translation order for a plus-strand gene is ascending across sorted parts
    assert plus.genomic_positions == [3, 4, 5, 9, 10, 11]
    minus = gene_model.Gene("g", "chr1", strand=-1, parts=[(3, 6), (9, 12)])
    # minus-strand translation order is reversed
    assert minus.genomic_positions == [11, 10, 9, 5, 4, 3]


def test_gene_add_part_appends():
    gene = gene_model.Gene("g", "chr1", strand=1)
    assert gene.genomic_start is None
    gene.add_part(3, 6)
    gene.add_part(9, 12)
    assert gene.parts == [(3, 6), (9, 12)]
    assert (gene.genomic_start, gene.genomic_end) == (3, 12)


def test_genemodel_derives_four_sequences_including_introns():
    # two exons [3,6) + [9,12) flanking an intron [6,9); coding = ATG|TTT, intron CCC
    contig = "AAA" + "ATG" + "CCC" + "TTT" + "AAA"
    model = gene_model.GeneModel("g", "c", strand=1, parts=[(3, 6), (9, 12)])
    model.finalize(contig)
    assert model.ref_coding == "ATGTTT"          # spliced CDS (no intron)
    assert model.protein == "MF"
    assert model.rna == "AUGUUU"                  # spliced CDS as RNA
    assert model.dna == "ATGCCCTTT"               # full span, intron included
    assert model.revcomp_dna == sequence.reverse_complement(model.dna)
    assert model.codon(1) == "ATG" and model.aa_at(1) == "M"


def test_genemodel_minus_strand_sequences_are_coding_oriented():
    contig = "AAA" + "ATG" + "CCC" + "TTT" + "AAA"
    model = gene_model.GeneModel("g", "c", strand=-1, parts=[(3, 6), (9, 12)])
    model.finalize(contig)
    # dna is the full span on the coding (minus) strand; revcomp is the plus strand
    assert model.dna == sequence.reverse_complement("ATGCCCTTT")
    assert model.revcomp_dna == "ATGCCCTTT"
    assert model.ref_coding == sequence.reverse_complement("ATGTTT")


# --------------------------------------------------------------------------- #
# theiagene.lib.parsers -- RawGene streams, matching, BED
# --------------------------------------------------------------------------- #

def test_match_identifiers_normalize_is_query_aware():
    quals = {"product": ["lanosterol 14-alpha demethylase"], "gene": ["ERG11"]}
    matched, identifiers = parsers.match_identifiers(
        quals, ["lanosterol.14-alpha.demethylase"], ["product", "gene"],
        exact_match=True, normalize=True,
    )
    # normalize path returns the query term that matched, spaces<->dots aware
    assert matched == "lanosterol.14-alpha.demethylase"
    assert "ERG11" in identifiers


def test_match_identifiers_raw_returns_matching_identifier():
    quals = {"product": ["FKS1"]}
    # substring (normalize=False): query 'FKS' matches identifier 'FKS1',
    # and the identifier itself is returned
    matched, _ = parsers.match_identifiers(
        quals, ["FKS"], ["product"], exact_match=False, normalize=False
    )
    assert matched == "FKS1"
    none_match, _ = parsers.match_identifiers(
        quals, ["ERG11"], ["product"], exact_match=False, normalize=False
    )
    assert none_match is None


def test_resolve_contig_prefers_membership_and_falls_back():
    assert parsers.resolve_contig(["a", "b"], {"b"}, True, "VCF") == "b"
    # not required -> first candidate even when absent
    assert parsers.resolve_contig(["a", "b"], {"z"}, False, "VCF") == "a"
    with pytest.raises(KeyError, match="a and b not in VCF"):
        parsers.resolve_contig(["a", "b"], {"z"}, True, "VCF")


def test_iter_gff_raw_coalesces_multi_exon_cds(tmp_path):
    path = tmp_path / "multi.gff"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tCDS\t1\t6\t.\t+\t0\tID=cds-A;product=geneA\n"
        "chr1\t.\tCDS\t11\t16\t.\t+\t0\tID=cds-A;product=geneA\n"
    )
    raws = list(parsers.iter_gff_raw(str(path), "CDS"))
    assert len(raws) == 1
    assert raws[0].parts == [(0, 6), (10, 16)]
    assert raws[0].qualifiers["product"] == ["geneA"]  # accumulated as a list
    assert raws[0].contig_seq is None                  # no FASTA supplied


def test_parse_bed_genes_coalesces_rows_by_name(tmp_path):
    bed = tmp_path / "q.bed"
    bed.write_text("chr1\t0\t6\tgeneA\nchr1\t10\t16\tgeneA\nchr1\t0\t3\tgeneB\n")
    genes = parsers.parse_bed_genes(str(bed), ["geneA", "geneB"], True, {"chr1"})
    by_id = {g.gene_id: g for g in genes}
    assert by_id["geneA"].parts == [(0, 6), (10, 16)]
    assert by_id["geneB"].parts == [(0, 3)]
