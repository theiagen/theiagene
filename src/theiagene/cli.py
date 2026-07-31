"""Unified ``theiagene`` command-line entrypoint.

Dispatches subcommands, each
of which registers its own arguments and runs its own pipeline."""

import sys
import argparse

from theiagene import __version__, gene_coverage, extract_variants
from theiagene.lib.logging_config import configure_logging


# subcommand name -> (module, help text)
_SUBCOMMANDS = (
    (
        "gene_coverage",
        gene_coverage,
        "quantify breadth/depth of coverage over query genes from a BAM",
    ),
    (
        "extract_variants",
        extract_variants,
        "extract a sub-VCF of variants within query genes from a VCF",
    ),
)


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level parser with one subparser per subcommand"""
    parser = argparse.ArgumentParser(
        prog="theiagene",
        description="Theiagen gene manipulation toolkit (theiagene)"
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    subparsers.required = True
    for name, module, help_text in _SUBCOMMANDS:
        sub = subparsers.add_parser(name, help=help_text, description=help_text)
        module.add_arguments(sub)
        sub.set_defaults(_handler=module.run_cli)
    return parser


def main(argv=None) -> int:
    """Parse arguments and dispatch to the selected subcommand"""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging()
    return args._handler(args)


if __name__ == "__main__":
    sys.exit(main())
