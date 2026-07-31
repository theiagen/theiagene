"""Extract the variants that fall within query-gene coordinates.

Given a VCF and query-gene coordinates (from a reference GFF or a BED file),
this command writes a sub-VCF containing only the variants that overlap the
``feature_type`` (e.g. CDS) segments of the query genes, annotating each kept
record with the overlapping query name(s) in a ``GENE`` INFO field."""

import sys
import logging
import argparse
from collections import defaultdict

import pysam

from theiagene.lib.feature import FeatureCol
from theiagene.lib.parsers import assimilate_gff, import_vcf
from theiagene.lib.query import (
    ordered_query_genes,
    extract_queries_from_bed,
    match_query,
    sanitize_info_value,
)
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)


def input_error_handling(args: argparse.Namespace) -> None:
    """Handle incompatible input arguments"""
    if not args.query_genes and not args.bedfile:
        raise ValueError("'query_genes' or 'bedfile' required")


def _iter_descendants(feature):
    """Yield every descendant of ``feature``, depth-first (children of children
    included), so a query unit's CDS are reachable no matter how deep the
    gene -> RNA -> CDS hierarchy runs."""
    for descendant in feature.descendants:
        yield descendant
        yield from _iter_descendants(descendant)


def _split_qualifiers(raw: str) -> list:
    """Split a comma-/space-delimited qualifier string into individual keys,
    dropping empty tokens."""
    return raw.replace(",", " ").split()


def feature_identifiers(feature, qualifiers) -> list:
    """Collect the candidate name strings a query term may match against.

    A gene can be named on the query feature itself, on its parent (the gene
    record usually carries ``gene``/``Name``) or on its CDS descendants (which
    usually carry ``product``); every such identifier is gathered so a query
    matches regardless of which record holds the name. The attribute keys to
    read are supplied by ``qualifiers`` and matched case-insensitively."""
    identifiers = []
    related = [feature]
    if feature.parent is not None:
        related.append(feature.parent)
    related.extend(_iter_descendants(feature))
    wanted = {qualifier.lower() for qualifier in qualifiers}
    for related_feature in related:
        if related_feature.fid:
            identifiers.append(related_feature.fid)
        for key, value in related_feature.attributes.items():
            if value and key.lower() in wanted:
                identifiers.append(value)
    return identifiers


def _feature_label(feature) -> str:
    """Return a stable human-readable name for a query feature, preferring an
    explicit gene name over the raw feature id."""
    for key in ("gene", "Name", "product"):
        value = feature.attributes.get(key)
        if value:
            return value
    return feature.fid


def gff_query_ranges(
    features: FeatureCol,
    query_list: list,
    group_by: str,
    feature_type: str,
    feature_qualifiers: list,
    exact_match: bool,
) -> dict:
    """Flatten the ``feature_type`` coordinates of the selected query genes to
    ``{<CONTIG>: [(START, END, LABEL), ...]}`` (0-based, half-open).

    The query units are the ``group_by`` features (e.g. each RNA). When
    ``query_list`` is non-empty a unit is kept only if one of its identifiers
    matches a query term, and the matched term becomes the annotation label;
    otherwise every unit is kept and labelled by its own gene name."""
    contig2ranges = defaultdict(list)
    for feature in features[group_by]:
        if query_list:
            label = match_query(
                query_list, feature_identifiers(feature, feature_qualifiers), exact_match
            )
            if label is None:
                continue
        else:
            label = _feature_label(feature)
        contig = feature.seqid
        # only the requested-type subfeatures beneath this query unit (e.g. CDS)
        subfeatures = FeatureCol(list(_iter_descendants(feature)), group=False)
        for subfeature in subfeatures[feature_type]:
            start, end = subfeature.start, subfeature.end
            if end <= start:
                raise ValueError(
                    f"Invalid region for query '{label}' on contig '{contig}': "
                    f"start ({start}) must be < end ({end})"
                )
            contig2ranges[contig].append((start, end, label))
    return contig2ranges


def bed_query_ranges(bedfile: str, query_set: set) -> dict:
    """Flatten a BED file to ``{<CONTIG>: [(START, END, NAME), ...]}`` (0-based,
    half-open), keeping only rows whose name is in ``query_set`` (an empty set
    keeps every row).

    A BED region is used directly as a query coordinate; the name column (col 4)
    is both the filter key and the annotation label."""
    contig2ranges = defaultdict(list)
    with open(bedfile) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            data = line.split()
            if not data:
                continue
            contig, start, end, name = data[0], int(data[1]), int(data[2]), data[3]
            if query_set and name not in query_set:
                continue
            if end <= start:
                raise ValueError(
                    f"Invalid region for query '{name}' on contig '{contig}': "
                    f"start ({start}) must be < end ({end})"
                )
            contig2ranges[contig].append((start, end, name))
    return contig2ranges


def collapse_to_single_contig(contig2ranges: dict, contig: str) -> dict:
    """Re-file every range under a single ``contig`` (the ambiguous-contig case,
    where the sample was mapped to a reference whose contig name differs from
    the coordinate source)."""
    collapsed = defaultdict(list)
    for ranges in contig2ranges.values():
        collapsed[contig].extend(ranges)
    return collapsed


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
        feature_qualifiers = _split_qualifiers(args.feature_qualifier)
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
