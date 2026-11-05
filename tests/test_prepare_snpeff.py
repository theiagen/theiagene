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
    with pytest.raises(ValueError, match="--translation_table"):
        ps.codon_table_lines("cauris", gff)


@pytest.mark.parametrize("input_table,gff_table", [(7, None), (None, 7)])
def test_codon_table_lines_rejects_unsupported_table(tmp_path, input_table, gff_table):
    gff = _gff(tmp_path, ("chr1", gff_table))
    with pytest.raises(ValueError, match="translation table 7"):
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


# --------------------------------------------------------------------------- #
# write_snpeff_config: template edge cases
# --------------------------------------------------------------------------- #

def _write_config(tmp_path, template_text):
    template = tmp_path / "template.config"
    template.write_text(template_text)
    config = tmp_path / "snpEff.config"
    ps.write_snpeff_config(
        str(template), str(config), str(tmp_path / "snpeff_data"), "cauris",
        "Candidozyma auris", ["cauris.codonTable : Alternative_Yeast_Nuclear\n"],
    )
    return config.read_text().splitlines()


def test_write_snpeff_config_handles_template_without_trailing_newline(tmp_path):
    # the genome registration must start its own line, not run onto the last one
    lines = _write_config(tmp_path, _TEMPLATE.rstrip("\n"))
    assert lines[-3:] == [
        "hg19.genome : Human",
        "cauris.genome : Candidozyma auris",
        "cauris.codonTable : Alternative_Yeast_Nuclear",
    ]


def test_write_snpeff_config_adds_data_dir_missing_from_template(tmp_path):
    # without a data.dir line SnpEff would fall back to ./data/ beside the
    # config instead of where genes.gff was written
    template = "".join(
        line for line in _TEMPLATE.splitlines(keepends=True) if not line.startswith("data.dir")
    )
    lines = _write_config(tmp_path, template)
    assert lines.count(f"data.dir = {tmp_path / 'snpeff_data'}/") == 1
    assert "cauris.genome : Candidozyma auris" in lines


# --------------------------------------------------------------------------- #
# database GFF pipe sanitization
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("line, expected, changed", [
    # NCBI gnl| identifiers lose their pipes in ID and Parent alike
    (
        "c1\tGenbank\tmRNA\t1\t9\t.\t-\t.\tID=rna-gnl|WGS:PGLS|mrna_T0;Parent=gene-1\n",
        "c1\tGenbank\tmRNA\t1\t9\t.\t-\t.\tID=rna-gnl_WGS:PGLS_mrna_T0;Parent=gene-1\n",
        True,
    ),
    (
        "c1\tGenbank\tCDS\t1\t9\t.\t-\t0\tID=cds-1;Parent=rna-gnl|WGS:PGLS|mrna_T0\n",
        "c1\tGenbank\tCDS\t1\t9\t.\t-\t0\tID=cds-1;Parent=rna-gnl_WGS:PGLS_mrna_T0\n",
        True,
    ),
    # comments, directives and pipe-free lines pass through untouched
    ("#comment|with|pipes\n", "#comment|with|pipes\n", False),
    ("c1\t.\tgene\t1\t9\t.\t+\t.\tID=gene-1\n", "c1\t.\tgene\t1\t9\t.\t+\t.\tID=gene-1\n", False),
])
def test_sanitize_gff_line(line, expected, changed):
    assert ps.sanitize_gff_line(line) == (expected, changed)


def test_main_strips_pipes_from_database_gff(tmp_path, caplog):
    reference = tmp_path / "reference.gff"
    reference.write_text(
        "##gff-version 3\n"
        "c1\tGenbank\tmRNA\t1\t300\t.\t+\t.\tID=rna-gnl|WGS:X|mrna_T0\n"
        "c1\tGenbank\tCDS\t1\t300\t.\t+\t0\tID=cds-1;Parent=rna-gnl|WGS:X|mrna_T0;transl_table=12\n"
    )
    template = tmp_path / "template.config"
    template.write_text(_TEMPLATE)
    data_dir = tmp_path / "snpeff_data"
    with caplog.at_level(logging.WARNING):
        ps.main([
            "--reference_gff", str(reference),
            "--template_config", str(template),
            "--data_dir", str(data_dir),
            "--genome_id", "cauris",
            "--organism", "Candidozyma auris",
            "--output", str(tmp_path / "snpEff.config"),
        ])
    gff = (data_dir / "cauris" / "genes.gff").read_text()
    assert "|" not in gff
    assert "Parent=rna-gnl_WGS:X_mrna_T0" in gff
    assert "2 GFF feature(s)" in caplog.text
