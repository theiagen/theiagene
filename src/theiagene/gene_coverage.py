"""Quantify breadth and depth of coverage over query genes.

Given a BAM and query-gene coordinates (from a reference GFF or a BED file),
this command reports the average depth, percent coverage, number of mapped
reads, and whether that read total meets ``--min_reads_mapped`` for each query
gene as JSON (``DEPTH_DICT.json``, ``COVERAGE_DICT.json``, ``READS_DICT.json``,
``READS_PASS_DICT.json``) and a readable ``COVERAGE_STATS.tsv``.

Query genes are selected from the GFF by matching ``--query_genes`` (or the BED
names) against each candidate feature's identifiers -- the ``--feature_qualifier``
attribute(s) (default ``product``) plus the feature/parent/descendant ids -- and
depth/coverage/reads are compiled from the matched unit's ``--feature_type``
(default ``CDS``) segments."""

import sys
import logging
import argparse
from collections import defaultdict

import pysam

from theiagene.lib.parsers import (
    write_json,
    import_bam,
    assimilate_gff,
)
from theiagene.lib.query import (
    ordered_query_genes,
    extract_queries_from_bed,
    split_qualifiers,
    gff_query_ranges,
    bed_query_ranges,
    collapse_to_single_contig,
    validate_region,
)
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)


def input_error_handling(args: argparse.Namespace) -> None:
    """Handle incompatible input arguments"""
    if not args.bedfile and not args.reference_gff:
        raise FileNotFoundError("'reference_gff' or 'bedfile' is required for coordinates")
    elif not args.query_genes and not args.bedfile:
        raise FileNotFoundError("'query_genes' or 'bedfile' required")


def union_intervals(intervals: list) -> list:
    """Merge overlapping or contiguous half-open ``(start, end)`` intervals into
    a minimal set of disjoint intervals, so a base shared by two ranges (e.g.
    overlapping isoform CDS) is represented once."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def count_mapped_reads(
    imported_bam: pysam.AlignmentFile,
    contig: str,
    start: int,
    end: int,
    min_mapping_quality: int = 0,
) -> set:
    """Identify the reads mapped over ``[start, end)`` on ``contig``.

    Returns a set of alignment identifiers rather than a count so the caller can
    pool ranges sharing a label (e.g. all CDS segments of one query gene) and
    count a read spanning several of them once. QC_fail/duplicate/secondary/
    supplementary alignments are excluded to match the reads ``count_coverage``
    feeds into depth and breadth, and ``min_mapping_quality`` filters reads by
    alignment confidence as ``min_quality`` filters bases by base quality."""
    read_ids = set()
    for read in imported_bam.fetch(contig, start, end):
        if (read.is_unmapped or read.is_secondary or read.is_supplementary
                or read.is_qcfail or read.is_duplicate):
            continue
        if read.mapping_quality < min_mapping_quality:
            continue
        # flag distinguishes the mates of a pair sharing one query name
        read_ids.add((read.query_name, read.flag, read.reference_start))
    return read_ids


def quantify_gene_coverage(
    imported_bam: pysam.AlignmentFile,
    contig2ranges: dict,
    min_depth: int = 1,
    min_quality: int = 0,
    expected_labels: list = None,
    min_mapping_quality: int = 0,
) -> tuple:
    """Quantify breadth/depth of coverage and mapped reads over query-gene coordinates.

    ``contig2ranges`` maps each contig to ``[(START, END, LABEL), ...]`` (0-based,
    half-open) as produced by :func:`gff_query_ranges`/:func:`bed_query_ranges`;
    every range sharing a ``LABEL`` (e.g. all CDS segments of one query gene) is
    pooled into a per-base union before counting, so a base covered by two
    overlapping ranges (e.g. isoform CDS) is counted once rather than twice. A
    read overlapping several of a label's ranges is likewise counted once.

    ``expected_labels`` (typically the resolved query list) seeds every output
    dict with 0 so a requested query that resolved to no coordinates is still
    reported as absent rather than silently dropped."""
    reference_names = set(imported_bam.references)
    # label -> contig -> [(start, end), ...], validated as we group
    label_ranges = defaultdict(lambda: defaultdict(list))

    for contig, ranges in contig2ranges.items():
        if contig not in reference_names:
            raise ValueError(f"Contig '{contig}' not found in BAM references")
        contig_len = imported_bam.get_reference_length(contig)
        for start, end, label in ranges:
            validate_region(start, end, label, contig, contig_len)
            label_ranges[label][contig].append((start, end))

    depth_dict = {label: 0 for label in expected_labels or []}
    coverage_dict = {label: 0 for label in expected_labels or []}
    reads_dict = {label: 0 for label in expected_labels or []}
    for label, contig_map in label_ranges.items():
        total_bases = 0
        summed_depth = 0
        covered_bases = 0
        # pooled across this label's ranges so a spanning read counts once
        read_ids = set()
        for contig, intervals in contig_map.items():
            # union overlapping ranges so shared bases are counted once
            for start, end in union_intervals(intervals):
                # excludes QC_fail/duplicate/secondary/supplementary reads
                coverage_data = imported_bam.count_coverage(
                    contig, start, end, quality_threshold=min_quality
                )
                read_ids |= count_mapped_reads(
                    imported_bam, contig, start, end, min_mapping_quality
                )
                for i, _ in enumerate(range(start, end)):
                    # calculate total depth across bases
                    base_depth = (
                        coverage_data[0][i]
                        + coverage_data[1][i]
                        + coverage_data[2][i]
                        + coverage_data[3][i]
                    )
                    summed_depth += base_depth
                    total_bases += 1
                    # base is considered covered if beyond minimum depth
                    if base_depth >= min_depth:
                        covered_bases += 1
        depth_dict[label] = summed_depth / total_bases
        # breadth is percent of covered bases exceeding min_depth
        coverage_dict[label] = 100 * (covered_bases / total_bases)
        # reads mapped anywhere across this label's pooled ranges
        reads_dict[label] = len(read_ids)

    return (
        dict(sorted(depth_dict.items())),
        dict(sorted(coverage_dict.items())),
        dict(sorted(reads_dict.items())),
    )


def flag_reads_mapped(reads_dict: dict, min_reads_mapped: int = 1) -> dict:
    """Flag each query gene as meeting ``min_reads_mapped`` mapped reads.

    The threshold is inclusive and applied to the reported totals, so the flag
    is a call layered onto the measurements rather than a filter of them -- a
    gene below the threshold keeps its measured depth, breadth, and reads."""
    return {
        query: reads >= min_reads_mapped for query, reads in reads_dict.items()
    }


def make_tsv(
    depth_dict: dict,
    coverage_dict: dict,
    reads_dict: dict,
    pass_dict: dict,
    ambiguous_contig: bool,
) -> str:
    """Make a readable TSV to convey depth, coverage, and mapped reads"""
    if ambiguous_contig:
        name = "query (WARNING: results may be inaccurate if sample is not mapped to reference used to generate BED file coordinates)"
    else:
        name = "query"
    tsv_str = f"#{name}\taverage_depth\tpercent_coverage\treads_mapped\treads_mapped_pass\n"
    for query, depth in depth_dict.items():
        tsv_str += (
            f"{query}\t{depth}\t{coverage_dict[query]}\t{reads_dict[query]}"
            f"\t{pass_dict[query]}\n"
        )
    return tsv_str


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the gene_coverage arguments on ``parser``"""
    parser.add_argument("--bam", required=True)
    parser.add_argument("--bedfile")
    parser.add_argument("--reference_gff")
    parser.add_argument("--query_genes", nargs="+")
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument(
        "--feature_qualifier",
        default="product",
        help="comma-/space-delimited attribute key(s), matched case-insensitively, "
        "used to collect query-gene name candidates",
    )
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument("--ambiguous_contig", action="store_true")
    parser.add_argument("--min_depth", type=int, default=1)
    parser.add_argument("--min_quality", type=int, default=0)
    parser.add_argument(
        "--min_mapping_quality",
        type=int,
        default=0,
        help="minimum alignment mapping quality for a read to count as mapped",
    )
    parser.add_argument(
        "--min_reads_mapped",
        type=int,
        default=1,
        help="minimum mapped reads for a query gene to pass (reported alongside "
        "the measurements, which are never filtered)",
    )
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Execute the gene_coverage pipeline and write output files"""
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
            args.feature_type,
            feature_qualifiers,
            args.exact_match,
        )
    else:
        contig2ranges = bed_query_ranges(args.bedfile, set(query_list))

    if not contig2ranges:
        logger.warning("No query-gene coordinates were resolved from the inputs")

    # import BAM and modify contig coordinates if needed
    imported_bam = import_bam(args.bam)

    if args.ambiguous_contig:
        # can't apply the ambiguous-contig approach with multiple contigs
        bam_contigs = imported_bam.references
        if len(bam_contigs) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple contigs in the BAM"
            )
        # re-file every query range under the single BAM contig
        contig2ranges = collapse_to_single_contig(contig2ranges, bam_contigs[0])

    # quantify statistics and write
    depth_dict, coverage_dict, reads_dict = quantify_gene_coverage(
        imported_bam, contig2ranges, args.min_depth, args.min_quality, query_list,
        args.min_mapping_quality
    )
    pass_dict = flag_reads_mapped(reads_dict, args.min_reads_mapped)
    write_json("DEPTH_DICT.json", depth_dict)
    write_json("COVERAGE_DICT.json", coverage_dict)
    write_json("READS_DICT.json", reads_dict)
    write_json("READS_PASS_DICT.json", pass_dict)

    tsv_str = make_tsv(
        depth_dict, coverage_dict, reads_dict, pass_dict,
        args.ambiguous_contig and args.bedfile
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
