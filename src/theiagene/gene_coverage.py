"""Quantify breadth and depth of coverage over query genes.

Given a BAM and query-gene coordinates (from a reference GBFF/GFF or a BED
file), this command reports the average depth and percent coverage of each
query gene as JSON (``DEPTH_DICT.json``, ``COVERAGE_DICT.json``) and a readable
``COVERAGE_STATS.tsv``."""

import sys
import logging
import argparse

import pysam

from theiagene.lib.query import extract_queries_from_bed
from theiagene.lib.gene_model import Gene
from theiagene.lib.parsers import (
    write_json,
    import_bam,
    iter_gbff_raw,
    iter_gff_raw,
    match_identifiers,
    resolve_contig,
    parse_bed_genes,
)
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)


def input_error_handling(args: argparse.Namespace) -> None:
    """Handle incompatible input arguments"""
    if not args.bedfile and not args.reference_gbff and not args.reference_gff:
        raise FileNotFoundError(
            "'reference_gbff', 'reference_gff', or 'bedfile' is required"
        )
    elif not args.query_genes and not args.bedfile:
        raise ValueError("'query_genes' or 'bedfile' required")


def _coverage_genes(
    raw_features,
    query_set: set,
    feature_qualifier: str,
    exact_match: bool,
    contig_names: set,
    require_contig: bool,
    feature_type: str = "CDS",
) -> list:
    """Turn a Gene stream into a list of matched Gene objects.

    Matching is the coverage flavour (raw exact/substring on the single feature
    qualifier); genes sharing a name on a contig accumulate their parts.  Only the
    ``feature_type`` coordinates are carried forward (filed under that type so
    :func:`quantify_gene_coverage` reads them back); a parsed gene carrying no
    part of that type (e.g. a non-coding gene when ``feature_type='CDS'``) is
    skipped."""
    qualifier = feature_qualifier.strip()
    genes = {}
    order = []
    for raw in raw_features:
        matched, _ = match_identifiers(
            raw.qualifiers, list(query_set), [qualifier],
            exact_match=exact_match, normalize=False,
        )
        segments = raw.segments(feature_type)
        if matched is None or not segments:
            continue
        contig = resolve_contig(raw.contig_candidates, contig_names, require_contig, "BAM")
        key = (contig, matched)
        gene = genes.get(key)
        if gene is None:
            gene = Gene(gene_id=matched, contig=contig, strand=raw.strand)
            genes[key] = gene
            order.append(key)
        else:
            logger.warning(f"{matched} recovered multiple times in {contig}")
        for start, end in segments:
            gene.add_part(start, end, feature=feature_type)
    return [genes[key] for key in order]


def quantify_gene_coverage(
    imported_bam: pysam.AlignmentFile,
    genes: list,
    feature_type: str = "CDS",
    min_depth: int = 1,
    min_quality: int = 0,
) -> tuple:
    """Quantify gene breadth and depth of coverage over a list of Gene objects.

    Breadth/depth are compiled from each gene's ``feature_type`` coordinate
    segments (parsed from its :attr:`~theiagene.lib.gene_model.Gene.parts`)."""
    depth_dict = {}
    coverage_dict = {}
    reference_names = set(imported_bam.references)

    for gene in genes:
        contig = gene.contig
        query = gene.gene_id
        if contig not in reference_names:
            raise ValueError(f"Contig '{contig}' not found in BAM references")
        contig_len = imported_bam.get_reference_length(contig)
        if query in depth_dict:
            logger.warning(
                f"{query} is present on multiple contigs and will be overwritten"
            )
        # check coverage data across range
        depths = []
        coverages = []
        for start, end in gene.segments(feature_type):
            start, end = int(start), int(end)
            if end <= start:
                raise ValueError(
                    f"Invalid region for query '{query}' on contig '{contig}': start ({start}) must be < end ({end})"
                )
            if start < 0:
                raise ValueError(
                    f"Invalid region for query '{query}' on contig '{contig}': start ({start}) must be >= 0"
                )
            if end > contig_len:
                raise ValueError(
                    f"Invalid region for query '{query}' on contig '{contig}': end ({end}) exceeds contig length ({contig_len})"
                )
            coverage_data = imported_bam.count_coverage(
                contig, start, end, quality_threshold=min_quality
            )
            for i, _ in enumerate(range(start, end)):
                # calculate total depth across bases
                total_depth = (
                    coverage_data[0][i]
                    + coverage_data[1][i]
                    + coverage_data[2][i]
                    + coverage_data[3][i]
                )
                # base is considered covered if beyond minimum depth
                coverages.append(total_depth >= min_depth)
                depths.append(total_depth)
        if not depths:
            raise ValueError(
                f"No positions evaluated for query '{query}' on contig '{contig}'"
            )
        depth_dict[query] = sum(depths) / len(depths)
        # breadth is percent of covered bases exceeding min_depth
        coverage_dict[query] = 100 * (sum(coverages) / len(coverages))

    return dict(sorted(depth_dict.items())), dict(sorted(coverage_dict.items()))


def make_tsv(depth_dict: dict, coverage_dict: dict, ambiguous_contig: bool) -> str:
    """Make a readable TSV to convey depth and coverage"""
    if ambiguous_contig:
        name = "query (WARNING: results may be inaccurate if sample is not mapped to reference used to generate BED file coordinates)"
    else:
        name = "query"
    tsv_str = f"#{name}\taverage_depth\tpercent_coverage\n"
    for query, depth in depth_dict.items():
        tsv_str += f"{query}\t{depth}\t{coverage_dict[query]}\n"
    return tsv_str.strip()


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the gene_coverage arguments on ``parser``"""
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gbff")
    parser.add_argument("--reference_gff")
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--feature_qualifier", default="product", help="feature qualifier to derive name search from")
    parser.add_argument(
        "--feature_type",
        default="CDS",
        help="feature type whose coordinates drive the coverage calculation "
        "(e.g. CDS, exon, gene); parsed from the reference GFF/GBFF (default: CDS)",
    )
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument("--ambiguous_contig", action="store_true")
    parser.add_argument("--min_depth", type=int, default=1)
    parser.add_argument("--min_quality", type=int, default=0)
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Execute the gene_coverage pipeline and write output files"""
    # error parsing
    input_error_handling(args)

    # import queries
    if args.query_genes:
        query_set = set()
        for queries in args.query_genes:
            query_set = query_set.union(q.strip() for q in queries.split(","))
    else:
        query_set = extract_queries_from_bed(args.bedfile)

    # import BAM and modify contig coordinates if needed
    imported_bam = import_bam(
        args.bam, args.ambiguous_contig
    )

    # collect query-gene coordinates as a list of Gene objects
    contig_names = set(imported_bam.references)
    require_contig = not args.ambiguous_contig
    genes = []
    if args.reference_gbff:
        genes = _coverage_genes(
            iter_gbff_raw(args.reference_gbff),
            query_set, args.feature_qualifier, args.exact_match,
            contig_names, require_contig, args.feature_type,
        )
    elif args.reference_gff:
        genes = _coverage_genes(
            iter_gff_raw(args.reference_gff),
            query_set, args.feature_qualifier, args.exact_match,
            contig_names, require_contig, args.feature_type,
        )
    if args.bedfile:
        genes += parse_bed_genes(
            args.bedfile, query_set, args.exact_match, contig_names,
            require_contig, feature_type=args.feature_type,
        )

    if args.ambiguous_contig:
        # reassign every gene to the single BAM contig
        contig = imported_bam.references[0]
        for gene in genes:
            gene.contig = contig

    # quantify statistics and write
    depth_dict, coverage_dict = quantify_gene_coverage(
        imported_bam, genes, args.feature_type, args.min_depth, args.min_quality
    )
    write_json("DEPTH_DICT.json", depth_dict)
    write_json("COVERAGE_DICT.json", coverage_dict)

    # add missing entries to TSV report
    missing_genes = query_set.difference(set(coverage_dict.keys()))
    for gene in missing_genes:
        # depth may be reported for those that have no breadth
        if gene not in depth_dict:
            depth_dict[gene] = 0
        coverage_dict[gene] = 0

    tsv_str = make_tsv(
        depth_dict, coverage_dict, args.ambiguous_contig and args.bedfile
    )
    with open("COVERAGE_STATS.tsv", "w") as out:
        out.write(tsv_str)

    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.gene_coverage``)"""
    parser = argparse.ArgumentParser(
        description="task_gene_coverage.wdl dependency script (github.com/theiagen/public_health_bioinformatics)"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
