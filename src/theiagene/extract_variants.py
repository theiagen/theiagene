"""Extract the variants that fall within query-gene coordinates.

Given a VCF and query-gene coordinates (from a reference GFF or a BED file),
this command writes a sub-VCF containing only the variants that overlap the
``feature_type`` (e.g. CDS) segments of the query genes, annotating each kept
record with the overlapping query name(s) in a ``GENE`` INFO field."""

import sys
import logging
import argparse

import pysam

from theiagene.lib.parsers import assimilate_gff, import_vcf
from theiagene.lib.query import (
    ordered_query_genes,
    extract_queries_from_bed,
    sanitize_info_value,
    split_qualifiers,
    gff_query_ranges,
    bed_query_ranges,
    collapse_to_single_contig,
)
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)


def input_error_handling(args: argparse.Namespace) -> None:
    """Handle incompatible input arguments"""
    if not args.bedfile and not args.reference_gff:
        raise FileNotFoundError("'reference_gff' or 'bedfile' is required for coordinates")
    elif not args.query_genes and not args.bedfile:
        raise ValueError("'query_genes' or 'bedfile' required")


def extract_variants(
    vcf_in: pysam.VariantFile, contig2ranges: dict, output_vcf: str
) -> int:
    """Write the variants overlapping the query coordinates to ``output_vcf``,
    annotating the overlapping query name(s) in a ``GENE`` INFO field.

    Returns the count of written records."""
    if "GENE" not in vcf_in.header.info:
        vcf_in.header.info.add(
            "GENE",
            ".",
            "String",
            "Query gene(s) whose extracted coordinate range overlaps this variant",
        )
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    written = 0
    for record in vcf_in:
        ranges = contig2ranges.get(record.contig)
        if not ranges:
            continue
        # pysam VariantRecord coordinates are 0-based, half-open (start, stop)
        overlapping = set(
            label
            for start, end, label in ranges
            if record.start < end and record.stop > start
        )
        if overlapping:
            clean = sorted(sanitize_info_value(label) for label in overlapping)
            record.info["GENE"] = clean
            vcf_out.write(record)
            written += 1

    vcf_out.close()
    vcf_in.close()
    logger.debug(f"Wrote {written} query-overlapping variant(s) to {output_vcf}")
    return written


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the extract_variants arguments on ``parser``"""
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gff")
    parser.add_argument("--group_by", default="RNA")
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument(
        "--feature_qualifier",
        default="Name,gene,product,locus_tag,Alias",
        help="comma-/space-delimited attribute key(s), matched case-insensitively, "
        "used to collect query-gene name candidates",
    )
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument("--ambiguous_contig", action="store_true")
    parser.add_argument("--output", default="EXTRACTED_VARIANTS.vcf")
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Extract variants associated with query genes and write the sub-VCF"""
    input_error_handling(args)

    # import queries: explicit --query_genes take precedence over BED names
    if args.query_genes:
        query_list = ordered_query_genes(args.query_genes)
    elif args.reference_gff:
        # BED supplies the query names when it is not the coordinate source
        query_list = sorted(extract_queries_from_bed(args.bedfile))
    else:
        # BED is itself the coordinate source; every row is a query
        query_list = []

    # build query coordinates from whichever source is authoritative
    if args.reference_gff:
        features = assimilate_gff(args.reference_gff)
        feature_qualifiers = split_qualifiers(args.feature_qualifier)
        contig2ranges = gff_query_ranges(
            features,
            query_list,
            args.group_by,
            args.feature_type,
            feature_qualifiers,
            args.exact_match,
        )
    else:
        contig2ranges = bed_query_ranges(args.bedfile, set(query_list))

    if not contig2ranges:
        logger.warning("No query-gene coordinates were resolved from the inputs")

    vcf_in = import_vcf(args.vcf)

    if args.ambiguous_contig:
        vcf_contigs = list(vcf_in.header.contigs)
        if len(vcf_contigs) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple "
                "contigs in the VCF"
            )
        contig_ranges = {tuple(r) for ranges in contig2ranges.values() for r in ranges}
        contig2ranges = collapse_to_single_contig(
            {0: contig_ranges}, vcf_contigs[0]
        )

    extract_variants(vcf_in, contig2ranges, args.output)

    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.extract_variants``)"""
    parser = argparse.ArgumentParser(
        description="create a sub-vcf extracting variants within query genes"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
