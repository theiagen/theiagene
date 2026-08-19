"""Tests for the unified ``theiagene`` command-line entrypoint wiring.

Only the parser wiring is exercised here: the gene_coverage subcommand's
end-to-end run path still has open design gaps (see the review notes), so it is
intentionally not driven through cli.main."""

import pytest

from theiagene import cli, gene_coverage


def test_build_parser_registers_gene_coverage_subcommand():
    parser = cli.build_parser()
    args = parser.parse_args(
        ["gene_coverage", "--bam", "x.bam", "--reference_gff", "r.gff"]
    )
    assert args._handler is gene_coverage.run_cli
    assert args.bam == "x.bam"
    assert args.reference_gff == "r.gff"


def test_gene_coverage_argument_defaults():
    parser = cli.build_parser()
    args = parser.parse_args(["gene_coverage", "--bam", "x.bam", "--reference_gff", "r.gff"])
    assert args.feature_type == "CDS"
    assert args.min_depth == 1
    assert args.min_base_quality == 0
    assert args.ambiguous_contig is False


def test_gene_coverage_accepts_deprecated_min_quality_alias():
    parser = cli.build_parser()
    args = parser.parse_args(
        ["gene_coverage", "--bam", "x.bam", "--reference_gff", "r.gff",
         "--min_quality", "20"]
    )
    assert args.min_base_quality == 20


def test_version_flag_exits_cleanly(capsys):
    with pytest.raises(SystemExit) as exc:
        cli.main(["--version"])
    assert exc.value.code == 0
    assert "theiagene" in capsys.readouterr().out


def test_missing_subcommand_errors():
    with pytest.raises(SystemExit) as exc:
        cli.main([])
    assert exc.value.code != 0
