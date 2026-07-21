"""VCF coordinate helpers shared by the theiagene commands.

``extract_vcf_genes`` filters a VCF to variants overlapping query-gene
coordinates; ``gene_coverage`` uses it to pre-extract a ``GENE_VARIANTS.vcf``
and ``variant_annotation`` consumes the resulting ``GENE`` INFO field.  Keeping
it here (rather than importing across command modules) gives both commands a
single, shared implementation of the coordinate model.  Both entry points take a
list of :class:`~theiagene.lib.gene_model.Gene` objects."""

import logging
from collections import defaultdict

import pysam

from theiagene.lib.query import sanitize_info_value


logger = logging.getLogger(__name__)


def flatten_coords_by_contig(genes, full_range: bool = False) -> dict:
    """Flatten Gene objects to {<CONTIG>: [(START, END, GENE_ID), ...]} for interval
    overlap testing.  If full_range is True, each gene is collapsed to a single
    (genomic_start, genomic_end) range spanning all of its parts; otherwise every
    part is emitted separately"""
    contig2ranges = defaultdict(list)
    for gene in genes:
        if full_range:
            # collapse all parts of a gene into a single spanning range
            contig2ranges[gene.contig].append(
                (gene.genomic_start, gene.genomic_end, gene.gene_id)
            )
        else:
            for start, end in gene.parts:
                contig2ranges[gene.contig].append((start, end, gene.gene_id))
    return contig2ranges


def extract_vcf_genes(vcffile: str, genes, output_vcf: str) -> int:
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

    # {<CONTIG>: [(START, END, GENE_ID), ...]} (0-based, half-open coordinates)
    contig2ranges = flatten_coords_by_contig(genes)

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
