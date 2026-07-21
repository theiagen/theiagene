"""VCF coordinate helpers shared by the theiagene commands.

``extract_vcf_genes`` filters a VCF to variants overlapping query-gene
coordinates; ``gene_coverage`` uses it to pre-extract a ``GENE_VARIANTS.vcf``
and ``variant_annotation`` consumes the resulting ``GENE`` INFO field.  Keeping
it here (rather than importing across command modules) gives both commands a
single, shared implementation of the coordinate model."""

import logging
from collections import defaultdict
from itertools import chain

import pysam

from theiagene.lib.query import sanitize_info_value


logger = logging.getLogger(__name__)


def flatten_coords_by_contig(contig2query2coords: dict, full_range: bool = False) -> dict:
    """Flatten to {<CONTIG>: [(START, END, QUERY), ...]} for interval overlap testing.
    If full_range is True, each query is collapsed to a single (min START, max END) range
    spanning all of its parts; otherwise every part is emitted separately"""
    contig2ranges = defaultdict(list)
    for contig, query2coords in contig2query2coords.items():
        for query, loc_parts in query2coords.items():
            if full_range:
                # collapse all parts of a query into a single spanning range
                all_coords = [int(coord) for coord in chain.from_iterable(loc_parts)]
                min_coord = min(all_coords)
                max_coord = max(all_coords)
                contig2ranges[contig].append((min_coord, max_coord, query))
            else:
                for coords in loc_parts:
                    contig2ranges[contig].append(
                        (int(coords[0]), int(coords[1]), query)
                    )
    return contig2ranges


def extract_vcf_genes(
    vcffile: str, contig2query2coords: dict, output_vcf: str
) -> int:
    """Filter a VCF to variants overlapping query gene coordinates, annotating the
    overlapping gene name(s) in a GENE INFO field. Returns the count of written records"""
    vcf_in = pysam.VariantFile(vcffile)
    # define the INFO field used to annotate the overlapping gene name(s)
    if "GENE" not in vcf_in.header.info:
        vcf_in.header.info.add(
            "GENE",
            ".",
            "String",
            "Query gene(s) whose extracted coordinate range overlaps this variant",
        )
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # {<CONTIG>: [(START, END, QUERY), ...]} (0-based, half-open coordinates)
    contig2ranges = flatten_coords_by_contig(contig2query2coords)

    written = 0
    for record in vcf_in:
        ranges = contig2ranges.get(record.contig)
        if not ranges:
            continue
        # pysam VariantRecord coordinates are 0-based, half-open (record.start, record.stop)
        genes = set(
            gene
            for start, end, gene in ranges
            if record.start < end and record.stop > start
        )
        if genes:
            # dedupe while preserving order and sanitize for the INFO field
            clean_genes = sorted([
                sanitize_info_value(gene) for gene in dict.fromkeys(genes)
            ])
            record.info["GENE"] = clean_genes
            vcf_out.write(record)
            written += 1

    vcf_out.close()
    vcf_in.close()
    logger.debug(f"Wrote {written} gene-overlapping variant(s) to {output_vcf}")
    return written
