"""Quantify breadth and depth of coverage over query genes.

Given a BAM and query-gene coordinates (from a reference GBFF/GFF or a BED
file), this command reports the average depth and percent coverage of each
query gene as JSON (``DEPTH_DICT.json``, ``COVERAGE_DICT.json``) and a readable
``COVERAGE_STATS.tsv``.  When a VCF is supplied, the gene-overlapping variants
are additionally extracted to ``GENE_VARIANTS.vcf``.

This is the ``gene_coverage`` subcommand of the ``theiagene`` entrypoint; it was
formerly the standalone ``gene_coverage.py`` script (a
``task_gene_coverage.wdl`` dependency in
github.com/theiagen/public_health_bioinformatics)."""

import sys
import logging
import argparse

import pysam

# shared helpers, re-exported so callers (and tests) can reach them here
from theiagene.lib.query import (  # noqa: F401
    exact_check,
    substring_check,
    extract_queries_from_bed,
)
from theiagene.lib.io_utils import write_json  # noqa: F401
from theiagene.lib.gene_model import Gene
from theiagene.lib.parsers import (
    import_bam,
    iter_gbff_raw,
    iter_gff_raw,
    match_identifiers,
    resolve_contig,
    parse_bed_genes,
    extract_vcf_genes,  # noqa: F401
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
) -> list:
    """Turn a Gene stream into a list of matched Gene objects.

    Matching is the coverage flavour (raw exact/substring on the single feature
    qualifier); genes sharing a name on a contig accumulate their parts.  A parsed
    gene carrying no CDS coordinates (e.g. a non-coding GFF gene) is skipped."""
    qualifier = feature_qualifier.strip()
    genes = {}
    order = []
    for raw in raw_features:
        matched, _ = match_identifiers(
            raw.qualifiers, list(query_set), [qualifier],
            exact_match=exact_match, normalize=False,
        )
        if matched is None or not raw.cds:
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
        for start, end in raw.cds:
            gene.add_part(start, end)
    return [genes[key] for key in order]


def quantify_gene_coverage(
    imported_bam: pysam.AlignmentFile,
    genes: list,
    min_depth: int = 1,
    min_quality: int = 0,
) -> tuple:
    """Quantify gene breadth and depth of coverage over a list of Gene objects"""
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
        for start, end in gene.cds:
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
    parser.add_argument("--vcf")
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gbff")
    parser.add_argument("--reference_gff")
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--feature_qualifier", default="product")
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
            contig_names, require_contig,
        )
    elif args.reference_gff:
        genes = _coverage_genes(
            iter_gff_raw(args.reference_gff),
            query_set, args.feature_qualifier, args.exact_match,
            contig_names, require_contig,
        )
    if args.bedfile:
        genes += parse_bed_genes(
            args.bedfile, query_set, args.exact_match, contig_names, require_contig
        )

    if args.ambiguous_contig:
        # reassign every gene to the single BAM contig
        contig = imported_bam.references[0]
        for gene in genes:
            gene.contig = contig

    # optionally extract gene-overlapping variants from a VCF into a single VCF;
    # the extraction routine lives in theiagene.lib.parsers (shared with the
    # variant_annotation command) and is re-exported into this module
    if args.vcf:
        extract_vcf_genes(args.vcf, genes, "GENE_VARIANTS.vcf")

    # quantify statistics and write
    depth_dict, coverage_dict = quantify_gene_coverage(
        imported_bam, genes, args.min_depth, args.min_quality
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
