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
from theiagene.lib.feature import FeatureCol
from theiagene.lib.parsers import (
    write_json,
    import_bam,
    assimilate_gff,
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


def quantify_gene_coverage(
    imported_bam: pysam.AlignmentFile,
    query_features: list,
    feature_type: str = "CDS",
    min_depth: int = 1,
    min_quality: int = 0,
) -> tuple:
    """Quantify gene breadth and depth of coverage over a list of query Features.

    Each Feature in ``query_features`` is a query unit (the ``group_by`` class,
    e.g. an RNA); breadth/depth are compiled from its ``feature_type`` (e.g. CDS)
    descendant segments. The contig is read from each Feature's ``seqid``."""
    depth_dict = {}
    coverage_dict = {}
    reference_names = set(imported_bam.references)

    for feature in query_features:
        contig = feature.seqid
        if contig not in reference_names:
            raise ValueError(f"Contig '{contig}' not found in BAM references")
        contig_len = imported_bam.get_reference_length(contig)
        query = feature.fid
        if query in depth_dict:
            logger.warning(
                f"{query} is present on multiple contigs and will be overwritten"
            )
        # check coverage data across range
        depths = []
        coverages = []
        # only the requested-type subfeatures beneath this query feature
        for subfeature in FeatureCol(feature.descendants, group=False)[feature_type]:
            start, end = subfeature.start, subfeature.end
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
    parser.add_argument("--reference_gff", required=True)
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--group_by", default="RNA")
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument("--exact_match", default=False)
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
    imported_bam = import_bam(args.bam)

    # convert GFF into a FeatureCol (parent/descendant hierarchy wired up)
    features = assimilate_gff(args.reference_gff)

    contig_names = set(imported_bam.references)
    if args.ambiguous_contig:
        # can't apply ambiguous contig approach if there are multiple contigs
        gff_contigs = {feature.seqid for feature in features}
        if len(contig_names) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple contigs in the BAM"
            )
        if len(gff_contigs) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple contigs in the GFF"
            )
        # reassign every feature to the single BAM contig
        only_contig = imported_bam.references[0]
        for feature in features:
            feature.seqid = only_contig

    # the group_by features (e.g. each RNA) are the query units
    query_features = features[args.group_by]

    # quantify statistics and write
    depth_dict, coverage_dict = quantify_gene_coverage(
        imported_bam, query_features, args.feature_type, args.min_depth, args.min_quality
    )
    write_json("DEPTH_DICT.json", depth_dict)
    write_json("COVERAGE_DICT.json", coverage_dict)

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
