"""Unit tests for theiagene.prepare_snpeff (SnpEff database GFF and config)."""

import logging

import pytest

from theiagene import prepare_snpeff as ps


_TEMPLATE = "\n".join(
    [
        "# SnpEff configuration file",
        "data.dir = ./data/",
        "codon.Standard : TTT/F",
        "hg19.genome : Human",
    ]
) + "\n"


def _gff(tmp_path, *tables, fasta=False):
    """A GFF with one CDS per (contig, transl_table) pair (table None omits it)"""
    rows = ["##gff-version 3"]
    for idx, (contig, table) in enumerate(tables):
        attributes = f"ID=cds-{idx}" + (f";transl_table={table}" if table else "")
        rows.append(f"{contig}\t.\tCDS\t1\t300\t.\t+\t0\t{attributes}")
    if fasta:
        rows += ["##FASTA", ">chr1", "ACGT"]
    path = tmp_path / "reference.gff"
    path.write_text("\n".join(rows) + "\n")
    return str(path)


# --------------------------------------------------------------------------- #
# codon_table_lines
# --------------------------------------------------------------------------- #

def test_codon_table_lines_assigns_each_contig_its_detected_table(tmp_path, caplog):
    gff = _gff(tmp_path, ("chr1", 12), ("chr2", 12), ("chrM", 3))
    with caplog.at_level(logging.DEBUG):
        lines = ps.codon_table_lines("cauris", gff)
    assert lines == [
        "cauris.chr1.codonTable : Alternative_Yeast_Nuclear\n",
        "cauris.chr2.codonTable : Alternative_Yeast_Nuclear\n",
        "cauris.chrM.codonTable : Yeast_Mitochondrial\n",
    ]
    assert "translation table 12 (Alternative_Yeast_Nuclear) for 2 contig(s): chr1, chr2" in caplog.text


def test_codon_table_lines_input_takes_precedence_genome_wide(tmp_path):
    # the GFF's conflicting tables are never consulted
    gff = _gff(tmp_path, ("chr1", 12), ("chr1", 3))
    assert ps.codon_table_lines("cauris", gff, 1) == ["cauris.codonTable : Standard\n"]


def test_codon_table_lines_warns_and_defaults_without_declarations(tmp_path, caplog):
    gff = _gff(tmp_path, ("chr1", None))
    with caplog.at_level(logging.WARNING):
        assert ps.codon_table_lines("cauris", gff) == []
    assert "no translation table detected" in caplog.text


def test_codon_table_lines_rejects_conflict_within_contig(tmp_path):
    gff = _gff(tmp_path, ("chr1", 12), ("chr1", 3))
    with pytest.raises(ValueError, match="must be input directly via --translation_table"):
        ps.codon_table_lines("cauris", gff)


@pytest.mark.parametrize("input_table,gff_table", [(7, None), (None, 7)])
def test_codon_table_lines_rejects_unsupported_table(tmp_path, input_table, gff_table):
    gff = _gff(tmp_path, ("chr1", gff_table))
    with pytest.raises(ValueError, match="translation table 7 is not supported"):
        ps.codon_table_lines("cauris", gff, input_table)


# --------------------------------------------------------------------------- #
# run_cli (end to end)
# --------------------------------------------------------------------------- #

def test_main_writes_database_gff_and_config(tmp_path):
    template = tmp_path / "template.config"
    template.write_text(_TEMPLATE)
    data_dir = tmp_path / "snpeff_data"
    config = tmp_path / "snpEff.config"
    ps.main([
        "--reference_gff", _gff(tmp_path, ("chr1", 12), fasta=True),
        "--template_config", str(template),
        "--data_dir", str(data_dir),
        "--genome_id", "cauris",
        "--organism", "Candidozyma auris",
        "--output", str(config),
    ])
    # the embedded FASTA section is dropped from the database GFF
    gff = (data_dir / "cauris" / "genes.gff").read_text()
    assert gff == "##gff-version 3\nchr1\t.\tCDS\t1\t300\t.\t+\t0\tID=cds-0;transl_table=12\n"
    assert config.read_text().splitlines() == [
        "# SnpEff configuration file",
        f"data.dir = {data_dir}/",
        "codon.Standard : TTT/F",
        "hg19.genome : Human",
        "cauris.genome : Candidozyma auris",
        "cauris.chr1.codonTable : Alternative_Yeast_Nuclear",
    ]


def test_main_writes_nothing_when_table_undeterminable(tmp_path):
    template = tmp_path / "template.config"
    template.write_text(_TEMPLATE)
    data_dir = tmp_path / "snpeff_data"
    with pytest.raises(ValueError):
        ps.main([
            "--reference_gff", _gff(tmp_path, ("chr1", 12), ("chr1", 3)),
            "--template_config", str(template),
            "--data_dir", str(data_dir),
            "--genome_id", "cauris",
            "--organism", "Candidozyma auris",
            "--output", str(tmp_path / "snpEff.config"),
        ])
    assert not data_dir.exists()
    assert not (tmp_path / "snpEff.config").exists()
