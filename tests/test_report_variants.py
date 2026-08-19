"""Unit tests for theiagene.report_variants (per-allele read-depth reporting)."""

import textwrap

import pytest

from theiagene import report_variants as rv
from theiagene.lib.parsers import assimilate_gff


# A minimal single-sample FreeBayes-style VCF exercising: a biallelic SNP, a
# biallelic indel (raw REF/ALT differ from VEP's normalized allele), a
# multiallelic site, a scientific-notation GQ (which import_vcf must scrub), and
# a gVCF reference block carrying no AD/RO/AO.
_VCF = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##contig=<ID=chr1>
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
    ##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Ref obs">
    ##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alt obs">
    ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Min depth">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
    chr1\t100\t.\tT\tC\t50\tPASS\tDP=562\tGT:GQ:DP:AD:RO:AO\t1:99:562:0,562:0:562
    chr1\t200\t.\tGA\tGAA\t40\tPASS\tDP=25\tGT:GQ:DP:AD:RO:AO\t1:3.5e+01:25:2,23:2:23
    chr1\t300\t.\tA\tG,T\t40\tPASS\tDP=30\tGT:GQ:DP:AD:RO:AO\t1:60:30:5,15,10:5:15,10
    chr1\t400\t.\tC\t<*>\t0\tPASS\tDP=44\tGT:GQ:DP:MIN_DP\t0:99:44:40
    """
)


@pytest.fixture
def vcf(tmp_path):
    path = tmp_path / "variants.vcf"
    path.write_text(_VCF)
    return str(path)


# --------------------------------------------------------------------------- #
# build_depth_index
# --------------------------------------------------------------------------- #

def test_build_depth_index_maps_raw_coordinates_to_allele_depths(vcf):
    index = rv.build_depth_index(vcf)
    assert index[("chr1", 100, "T")] == (0, {"C": 562})
    # scientific-notation GQ was scrubbed rather than raising
    assert index[("chr1", 200, "GA")] == (2, {"GAA": 23})
    # multiallelic keeps a depth per alt
    assert index[("chr1", 300, "A")] == (5, {"G": 15, "T": 10})


def test_build_depth_index_skips_records_without_allele_depths(vcf):
    # the gVCF reference block (no AD/RO/AO) is absent
    assert ("chr1", 400, "C") not in rv.build_depth_index(vcf)


# --------------------------------------------------------------------------- #
# _depth_suffix
# --------------------------------------------------------------------------- #

def test_depth_suffix_renders_ref_then_alt(vcf):
    index = rv.build_depth_index(vcf)
    row = {"Uploaded_variation": "chr1_100_T/C", "Allele": "C"}
    assert rv._depth_suffix(row, index) == "T:0 C:562"


def test_depth_suffix_falls_back_to_sole_alt_when_allele_is_normalized(vcf):
    # VEP normalizes the indel's Allele to "A"; the sole raw alt is still chosen
    index = rv.build_depth_index(vcf)
    row = {"Uploaded_variation": "chr1_200_GA/GAA", "Allele": "A"}
    assert rv._depth_suffix(row, index) == "GA:2 GAA:23"


def test_depth_suffix_selects_matching_alt_of_multiallelic(vcf):
    index = rv.build_depth_index(vcf)
    row = {"Uploaded_variation": "chr1_300_A/G,T", "Allele": "T"}
    assert rv._depth_suffix(row, index) == "A:5 T:10"


def test_depth_suffix_none_when_multiallelic_allele_unresolved(vcf):
    index = rv.build_depth_index(vcf)
    row = {"Uploaded_variation": "chr1_300_A/G,T", "Allele": "-"}
    assert rv._depth_suffix(row, index) is None


def test_depth_suffix_none_when_variant_absent(vcf):
    index = rv.build_depth_index(vcf)
    row = {"Uploaded_variation": "chr1_999_A/C", "Allele": "C"}
    assert rv._depth_suffix(row, index) is None


def test_depth_suffix_none_on_unparseable_upload(vcf):
    index = rv.build_depth_index(vcf)
    assert rv._depth_suffix({"Uploaded_variation": "garbage"}, index) is None


# --------------------------------------------------------------------------- #
# report_line
# --------------------------------------------------------------------------- #

def test_report_line_appends_depths_after_semicolon(vcf):
    index = rv.build_depth_index(vcf)
    row = {
        "Uploaded_variation": "chr1_100_T/C",
        "Allele": "C",
        "Consequence": "missense_variant",
        "HGVSc": "rna-x:c.428A>G",
        "HGVSp": "prot-x:p.(Lys143Arg)",
    }
    line = rv.report_line(row, "lanosterol 14-alpha demethylase", index)
    assert line == (
        'lanosterol.14-alpha.demethylase: "lanosterol 14-alpha demethylase" '
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)"
    )


def test_report_line_leads_with_query_label(vcf):
    index = rv.build_depth_index(vcf)
    row = {
        "Uploaded_variation": "chr1_100_T/C",
        "Allele": "C",
        "Consequence": "missense_variant",
        "HGVSc": "rna-x:c.428A>G",
        "HGVSp": "prot-x:p.(Lys143Arg)",
    }
    line = rv.report_line(row, "lanosterol 14-alpha demethylase", index, "ERG11")
    assert line == (
        'ERG11: "lanosterol 14-alpha demethylase" '
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)"
    )


def test_report_line_fallback_label_drops_commas_from_product():
    # the report is joined with ',', so a comma may not survive into the label
    row = {"Consequence": "missense_variant", "HGVSc": "rna-x:c.1A>C", "HGVSp": "-"}
    line = rv.report_line(row, "1,3-beta-glucan synthase component FKS1")
    assert line == (
        '1.3-beta-glucan.synthase.component.FKS1: '
        '"1,3-beta-glucan synthase component FKS1" (missense_variant c.1A>C)'
    )


def test_report_line_without_index_is_unchanged():
    row = {
        "Consequence": "missense_variant",
        "HGVSc": "rna-x:c.428A>G",
        "HGVSp": "-",
    }
    line = rv.report_line(row, "product name")
    assert line == 'product.name: "product name" (missense_variant c.428A>G)'


def test_report_line_expands_percent_encoded_synonymous_change():
    row = {
        "Consequence": "synonymous_variant",
        "HGVSc": "rna-x:c.492C>T",
        # VEP percent-encodes the '=' of a synonymous protein change
        "HGVSp": "prot-x:p.Asp164%3D",
    }
    line = rv.report_line(row, "product")
    assert line == 'product: "product" (synonymous_variant c.492C>T p.Asp164Asp)'


@pytest.mark.parametrize(
    "suffix,expected",
    [
        ("p.Asp164=", "p.Asp164Asp"),
        ("p.D164=", "p.D164D"),
        ("p.Ter330=", "p.Ter330Ter"),
        # no single reference residue to repeat -- left alone
        ("p.=", "p.="),
        ("p.Asp164_Leu166=", "p.Asp164_Leu166="),
        # not a synonymous change
        ("p.Lys143Arg", "p.Lys143Arg"),
        ("c.428A>G", "c.428A>G"),
    ],
)
def test_expand_synonymous(suffix, expected):
    assert rv._expand_synonymous(suffix) == expected


def test_report_line_omits_depths_when_variant_absent(vcf):
    index = rv.build_depth_index(vcf)
    row = {
        "Uploaded_variation": "chr1_999_A/C",
        "Allele": "C",
        "Consequence": "missense_variant",
        "HGVSc": "rna-x:c.1A>C",
        "HGVSp": "-",
    }
    line = rv.report_line(row, "product", index)
    assert line == 'product: "product" (missense_variant c.1A>C)'


# --------------------------------------------------------------------------- #
# report_variants (end to end, through the GFF feature hierarchy)
# --------------------------------------------------------------------------- #

# gene -> mRNA -> CDS, the hierarchy a VEP 'Feature' (the mRNA) is resolved
# through; the product's ',' is percent-encoded as GFF3 column 9 requires
_GFF = "\n".join(
    [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t1000\t.\t+\t.\tID=gene-FKS1;Name=FKS1",
        "chr1\t.\tmRNA\t1\t1000\t.\t+\t.\tID=rna-x;Parent=gene-FKS1",
        "chr1\t.\tCDS\t1\t1000\t.\t+\t0\tID=cds-x;Parent=rna-x;"
        "product=1%2C3-beta-glucan synthase component FKS1",
    ]
) + "\n"

_VEP_TSV = "\n".join(
    [
        "## VEP run statistics",
        "#Uploaded_variation\tLocation\tAllele\tConsequence\tFeature\tHGVSc\tHGVSp",
        "chr1_100_T/C\tchr1:100\tC\tsynonymous_variant\trna-x\t"
        "rna-x:c.492C>T\tprot-x:p.Asp164%3D",
    ]
) + "\n"


@pytest.fixture
def annotations(tmp_path):
    gff = tmp_path / "reference.gff"
    gff.write_text(_GFF)
    tsv = tmp_path / "annotations.tsv"
    tsv.write_text(_VEP_TSV)
    return str(gff), str(tsv)


def test_report_variants_labels_lines_by_query_gene(annotations):
    gff, tsv = annotations
    features = assimilate_gff(gff)
    lines = rv.report_variants(
        tsv, features, set(), "CDS", ["product"], query_list=["FKS1"]
    )
    assert lines == [
        'FKS1: "1,3-beta-glucan synthase component FKS1" '
        "(synonymous_variant c.492C>T p.Asp164Asp)"
    ]


def test_report_variants_falls_back_to_product_label(annotations):
    gff, tsv = annotations
    features = assimilate_gff(gff)
    # no query matches this feature, so the product names the line instead
    lines = rv.report_variants(
        tsv, features, set(), "CDS", ["product"], query_list=["ERG11"]
    )
    assert lines == [
        '1.3-beta-glucan.synthase.component.FKS1: '
        '"1,3-beta-glucan synthase component FKS1" '
        "(synonymous_variant c.492C>T p.Asp164Asp)"
    ]


def test_report_variants_drops_row_whose_hgvs_columns_are_absent(tmp_path):
    # VEP run without --hgvs emits no HGVSc/HGVSp columns at all; such a row
    # describes no variant, so it must be dropped rather than reported as a
    # bare consequence
    gff = tmp_path / "reference.gff"
    gff.write_text(_GFF)
    tsv = tmp_path / "no_hgvs.tsv"
    tsv.write_text(
        "#Uploaded_variation\tLocation\tAllele\tConsequence\tFeature\n"
        "chr1_100_T/C\tchr1:100\tC\tmissense_variant\trna-x\n"
    )
    features = assimilate_gff(str(gff))
    assert rv.report_variants(str(tsv), features, set(), "CDS", ["product"]) == []


def test_report_variants_drops_row_truncated_before_its_hgvs_columns(tmp_path):
    # a short data row leaves the trailing columns missing from the zipped dict,
    # which must read as undefined rather than as a value
    gff = tmp_path / "reference.gff"
    gff.write_text(_GFF)
    tsv = tmp_path / "ragged.tsv"
    tsv.write_text(
        "#Uploaded_variation\tAllele\tConsequence\tFeature\tHGVSc\tHGVSp\n"
        "chr1_100_T/C\tC\tmissense_variant\trna-x\n"
    )
    features = assimilate_gff(str(gff))
    assert rv.report_variants(str(tsv), features, set(), "CDS", ["product"]) == []


def test_report_variants_exact_match_rejects_substring_query(annotations):
    gff, tsv = annotations
    features = assimilate_gff(gff)
    # 'FKS1' is only a substring of the product/gene id, so --exact_match drops it
    lines = rv.report_variants(
        tsv, features, set(), "CDS", ["product"], query_list=["FKS1"], exact_match=True
    )
    assert lines[0].startswith("1.3-beta-glucan.synthase.component.FKS1: ")
