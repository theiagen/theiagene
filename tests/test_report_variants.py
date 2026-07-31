"""Unit tests for theiagene.report_variants (per-allele read-depth reporting)."""

import textwrap

import pytest

from theiagene import report_variants as rv


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
        "lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase "
        "(missense_variant c.428A>G p.Lys143Arg; T:0 C:562)"
    )


def test_report_line_without_index_is_unchanged():
    row = {
        "Consequence": "missense_variant",
        "HGVSc": "rna-x:c.428A>G",
        "HGVSp": "-",
    }
    line = rv.report_line(row, "product name")
    assert line == "product.name: product name (missense_variant c.428A>G)"


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
    assert line == "product: product (missense_variant c.1A>C)"
