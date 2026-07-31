"""Unit tests for theiagene.lib.parsers (JSON, GFF3 readers, BAM import)."""

import json

import pytest

from theiagene.lib import parsers
from theiagene.lib.feature import Feature


# --------------------------------------------------------------------------- #
# write_json
# --------------------------------------------------------------------------- #

def test_write_json_writes_data(tmp_path):
    out = tmp_path / "d.json"
    parsers.write_json(str(out), {"geneA": 5.0})
    assert json.loads(out.read_text()) == {"geneA": 5.0}


def test_write_json_spoofs_cromwell_on_empty(tmp_path):
    # an empty dict must still yield valid JSON for the WDL runner
    out = tmp_path / "empty.json"
    parsers.write_json(str(out), {})
    assert out.read_text() == '{"": 0}'


# --------------------------------------------------------------------------- #
# GFF3 attribute parsing / formatting
# --------------------------------------------------------------------------- #

def test_parse_gff_attributes_splits_and_percent_decodes():
    parsed = parsers.parse_gff_attributes(
        "ID=cds-1;product=lanosterol%2014-alpha;gene=ERG11;"
    )
    assert parsed == {
        "ID": "cds-1",
        "product": "lanosterol 14-alpha",  # %20 decoded to a space
        "gene": "ERG11",
    }


def test_parse_gff_attributes_empty_column_is_empty_dict():
    assert parsers.parse_gff_attributes("") == {}
    assert parsers.parse_gff_attributes(";") == {}


def test_parse_gff_attributes_raises_on_malformed_field():
    with pytest.raises(ValueError, match="unexpected attributes field"):
        parsers.parse_gff_attributes("novalue;ID=x")


def test_parse_gff_attributes_keeps_value_containing_delimiter():
    # only the first '=' splits key from value
    assert parsers.parse_gff_attributes("note=a=b") == {"note": "a=b"}


def test_format_gff_attributes_is_inverse_of_parse():
    assert parsers.format_gff_attributes({}) == "."
    assert parsers.format_gff_attributes({"ID": "x", "gene": "ERG11"}) == "ID=x;gene=ERG11"


# --------------------------------------------------------------------------- #
# iter_gff_features -> Feature stream
# --------------------------------------------------------------------------- #

def test_iter_gff_features_yields_features_with_converted_coordinates(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a;Parent=p;product=alpha\n"
        "chr2\t.\tCDS\t6\t17\t.\t-\t0\tID=b;product=beta\n"
    )
    feats = list(parsers.iter_gff_features(str(path)))
    assert len(feats) == 2
    first = feats[0]
    # GFF 1-based inclusive [10, 33] -> 0-based half-open [9, 33)
    assert (first.seqid, first.start, first.end, first.strand) == ("chr1", 9, 33, 1)
    # ID / Parent are lifted onto fid / pid and removed from the attribute dict
    assert first.fid == "a" and first.pid == "p"
    assert first.attributes == {"product": "alpha"}
    assert feats[1].strand == -1 and feats[1].pid is None


def test_iter_gff_features_skips_comments_blanks_and_fasta(tmp_path):
    path = tmp_path / "f.gff"
    path.write_text(
        "##gff-version 3\n"
        "\n"
        "# a comment\n"
        "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a\n"
        "##FASTA\n"
        "chr1\t.\tCDS\t1\t3\t.\t+\t0\tID=ignored\n"  # after ##FASTA: ignored
    )
    feats = list(parsers.iter_gff_features(str(path)))
    assert [f.fid for f in feats] == ["a"]


def test_iter_gff_features_reads_gzipped_gff(tmp_path):
    import gzip

    path = tmp_path / "f.gff.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "##gff-version 3\n"
            "# a comment\n"
            "chr1\t.\tCDS\t10\t33\t.\t+\t0\tID=a\n"
            "chr2\t.\tCDS\t6\t17\t.\t-\t0\tID=b\n"
        )
    # a .gz suffix enters the same parsing loop as a plain-text GFF
    feats = list(parsers.iter_gff_features(str(path)))
    assert [f.fid for f in feats] == ["a", "b"]


def test_iter_gff_features_raises_on_wrong_field_count(tmp_path):
    path = tmp_path / "bad.gff"
    path.write_text("chr1\t.\tCDS\t10\t33\t.\t+\n")  # 7 columns, not 9
    with pytest.raises(ValueError, match="incorrectly formatted GFF"):
        list(parsers.iter_gff_features(str(path)))


def test_iter_gff_features_raises_when_id_missing(tmp_path):
    path = tmp_path / "noid.gff"
    path.write_text("chr1\t.\tCDS\t10\t33\t.\t+\t0\tproduct=alpha\n")
    with pytest.raises(ValueError, match="found in record attributes"):
        list(parsers.iter_gff_features(str(path)))


# --------------------------------------------------------------------------- #
# assimilate_gff -> grouped hierarchy
# --------------------------------------------------------------------------- #

def test_assimilate_gff_groups_cds_under_its_rna_and_gene(tmp_path):
    path = tmp_path / "hier.gff"
    path.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t30\t.\t+\t.\tID=gene-A;gene=ERG11\n"
        "chr1\t.\tmRNA\t1\t30\t.\t+\t.\tID=rna-A;Parent=gene-A;product=geneA\n"
        "chr1\t.\tCDS\t5\t25\t.\t+\t0\tID=cds-A;Parent=rna-A;product=geneA\n"
    )
    grouped = parsers.assimilate_gff(str(path))
    assert set(grouped) == {"chr1"}
    by_id = {f.fid: f for f in grouped["chr1"]}
    rna = by_id["rna-A"]
    assert rna.parent.fid == "gene-A"
    # the CDS is a descendant of the mRNA, 1-based [5, 25] -> 0-based [4, 25)
    cds = rna.descendants[0]
    assert (cds.fid, cds.start, cds.end) == ("cds-A", 4, 25)


# --------------------------------------------------------------------------- #
# import_bam
# --------------------------------------------------------------------------- #

def test_import_bam_opens_and_indexes(make_bam):
    bam_path = make_bam(contig="contig1", contig_len=100)
    imported = parsers.import_bam(bam_path)
    assert "contig1" in set(imported.references)
    assert imported.get_reference_length("contig1") == 100
    assert imported.has_index()
