"""Prepare a reference GFF and config for building a custom SnpEff database.

SnpEff builds a database from ``<data_dir>/<genome_id>/genes.gff`` alongside the
reference sequences, under a config that registers the genome. This command:

- writes the GFF's annotation section to ``<data_dir>/<genome_id>/genes.gff``,
  dropping any embedded ``##FASTA`` section (SnpEff reads the sequences from
  the reference FASTA instead)
- writes a config derived from ``--template_config`` (retaining its codon tables
  and settings) that points ``data.dir`` at ``--data_dir``, registers the genome
  under ``--organism``, and assigns its codon table(s)

The codon table is ``--translation_table`` genome-wide when given; otherwise each
contig takes the ``transl_table`` its GFF features declare (see
:func:`theiagene.lib.parsers.gff_translation_tables`), and contigs declaring
none fall back to SnpEff's default (Standard) table."""

import os
import re
import sys
import logging
import argparse
from collections import defaultdict

from theiagene.lib.parsers import gff_translation_tables, iter_gff_lines
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)

# NCBI translation table number -> SnpEff codon table name
CODON_TABLES = {
    1: "Standard",
    2: "Vertebrate_Mitochondrial",
    3: "Yeast_Mitochondrial",
    4: "Mold_Mitochondrial",
    5: "Invertebrate_Mitochondrial",
    6: "Ciliate_Nuclear",
    9: "Echinoderm_Mitochondrial",
    10: "Euplotid_Nuclear",
    11: "Bacterial_and_Plant_Plastid",
    12: "Alternative_Yeast_Nuclear",
    13: "Ascidian_Mitochondrial",
    14: "Alternative_Flatworm_Mitochondrial",
    15: "Blepharisma_Macronuclear",
    16: "Chlorophycean_Mitochondrial",
    21: "Trematode_Mitochondrial",
    22: "Scenedesmus_obliquus_Mitochondrial",
    23: "Thraustochytrium_Mitochondrial",
}


def codon_table(table: int) -> str:
    """Return the SnpEff codon table name for an NCBI translation table number"""
    if table not in CODON_TABLES:
        raise ValueError(
            f"translation table {table} is not supported by SnpEff; supported "
            f"tables: {', '.join(map(str, CODON_TABLES))}"
        )
    return CODON_TABLES[table]


def codon_table_lines(genome_id: str, reference_gff: str, translation_table: int = None) -> list:
    """Return the config lines assigning ``genome_id`` its codon table(s).

    SnpEff reads a genome-wide codon table from ``<genome>.codonTable`` and a
    per-contig one from ``<genome>.<contig>.codonTable``. ``translation_table``
    takes precedence as a genome-wide table; otherwise each contig is assigned
    the table its GFF features declare, and no lines are returned when none
    declare one."""
    if translation_table is not None:
        name = codon_table(translation_table)
        logger.debug(f"using input translation table {translation_table} ({name}) for all contigs")
        return [f"{genome_id}.codonTable : {name}\n"]

    try:
        detected = gff_translation_tables(reference_gff)
    except ValueError as error:
        raise ValueError(
            f"{error}; the translation table cannot be determined and must be "
            "input directly via --translation_table"
        ) from None
    if not detected:
        logger.warning(
            "no translation table detected in the GFF; using the default (Standard) "
            "translation table"
        )
        return []

    contigs_by_table = defaultdict(list)
    lines = []
    for contig, table in detected.items():
        contigs_by_table[table].append(contig)
        lines.append(f"{genome_id}.{contig}.codonTable : {codon_table(table)}\n")
    for table, contigs in contigs_by_table.items():
        logger.debug(
            f"using GFF-detected translation table {table} ({CODON_TABLES[table]}) "
            f"for {len(contigs)} contig(s): {', '.join(contigs)}"
        )
    return lines


def write_snpeff_config(template_config: str, output: str, data_dir: str, genome_id: str,
                        organism: str, codon_lines: list) -> None:
    """Copy ``template_config`` to ``output``, pointing ``data.dir`` at
    ``data_dir`` and appending the genome registration and ``codon_lines``"""
    with open(template_config) as template, open(output, "w") as config:
        for line in template:
            if re.match(r"data\.dir\s*=", line):
                line = f"data.dir = {os.path.join(os.path.abspath(data_dir), '')}\n"
            config.write(line)
        config.write(f"{genome_id}.genome : {organism}\n")
        config.writelines(codon_lines)


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the prepare_snpeff arguments on ``parser``"""
    parser.add_argument("--reference_gff", required=True)
    parser.add_argument(
        "--template_config",
        required=True,
        help="SnpEff config the output config is derived from (e.g. the one "
        "installed alongside snpEff.jar)",
    )
    parser.add_argument(
        "--data_dir",
        required=True,
        help="SnpEff database directory; the GFF is written to "
        "<data_dir>/<genome_id>/genes.gff",
    )
    parser.add_argument(
        "--genome_id",
        required=True,
        help="SnpEff genome identifier (a config key, so restrict it to [A-Za-z0-9._-])",
    )
    parser.add_argument(
        "--organism",
        required=True,
        help="organism name the genome is registered under",
    )
    parser.add_argument(
        "--translation_table",
        type=int,
        help="NCBI translation table number applied genome-wide; detected per "
        "contig from the GFF transl_table attribute when unset",
    )
    parser.add_argument(
        "--output",
        default="snpEff.config",
        help="write the SnpEff config here",
    )
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Write the SnpEff database GFF and config"""
    codon_lines = codon_table_lines(args.genome_id, args.reference_gff, args.translation_table)

    genome_dir = os.path.join(args.data_dir, args.genome_id)
    os.makedirs(genome_dir, exist_ok=True)
    gff = os.path.join(genome_dir, "genes.gff")
    with open(gff, "w") as out:
        out.writelines(iter_gff_lines(args.reference_gff))
    logger.debug(f"Wrote the GFF annotation section to {gff}")

    write_snpeff_config(
        args.template_config, args.output, args.data_dir, args.genome_id,
        args.organism, codon_lines,
    )
    logger.debug(f"Wrote the SnpEff config to {args.output}")
    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.prepare_snpeff``)"""
    parser = argparse.ArgumentParser(
        description="Prepares a reference GFF and config for building a SnpEff database"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
