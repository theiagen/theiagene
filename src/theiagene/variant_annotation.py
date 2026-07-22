"""Annotate the protein-level consequences of variants that overlap query genes.

This is the ``variant_annotation`` subcommand of the ``theiagene`` entrypoint
(formerly the standalone ``variant_annotation.py`` script). It reports the
effect of every variant on the translated protein, given a VCF, a
reference GBFF/GFF + FA that supplies the gene coordinates and coding sequence, strand, product name
and translation table for each gene, and a set of query genes. 

Query genes are matched to each variant by interval overlap against the coding
models built from the reference. The
query-gene-overlapping variants are additionally written to a VCF
(``GENE_VARIANTS.vcf`` by default; ``--gene_vcf``) with a ``GENE`` INFO field
naming the overlapping gene(s).

For each variant it determines whether a substitution is

* ``synonymous_variant`` (silent, same amino acid),
* ``missense_variant``   (different amino acid),
* ``stop_gained``        (nonsense, a new stop codon),
* ``stop_lost``/``start_lost`` (edge substitutions),

and whether an indel is an ``inframe_deletion``, ``inframe_insertion`` or a
``frameshift_variant``.  A frameshift invalidates every downstream codon, so by
default all variants that fall 3' of a frameshift within the same gene are
disregarded; pass ``--annotate_downstream_of_frameshift`` (or
``suppress_downstream_frameshift=False``) to annotate them anyway.

The report is emitted as a single string, one entry per surviving variant::

    <gene_id>: <product> (<so_term> <hgvs_c> <hgvs_p>; <REF>:<depth> <ALT>:<depth>)

with entries joined by commas.  When no variants are recovered a single sentence
naming the queried genes is emitted instead.
"""

import sys
import logging
import argparse
from collections import defaultdict

import pysam

# shared helpers, re-exported so callers (and tests) can reach them here
from theiagene.lib.sequence import is_nucleotide_allele
from theiagene.lib.query import ordered_query_genes
from theiagene.lib.gene_model import GeneModel  # noqa: F401
from theiagene.lib.parsers import (  # noqa: F401
    build_gene_models_gbff,
    build_gene_models_gff,
    extract_vcf_genes,
)
from theiagene.lib.logging_config import configure_logging
from theiagene.lib.variant import Variant


logger = logging.getLogger(__name__)


def normalize_indel(pos0: int, ref: str, alt: str) -> tuple:
    """Trim shared prefix/suffix from a REF/ALT pair.

    Returns (changed_pos0, ref_segment, alt_segment) where ``changed_pos0`` is
    the 0-based genomic position of the first base of ``ref_segment``.  For a
    pure insertion ``ref_segment`` is empty and the insertion falls between
    ``changed_pos0 - 1`` and ``changed_pos0``; for a pure deletion
    ``alt_segment`` is empty."""
    ref = ref.upper()
    alt = alt.upper()
    # trim shared prefix
    prefix = 0
    while prefix < len(ref) and prefix < len(alt) and ref[prefix] == alt[prefix]:
        prefix += 1
    # trim shared suffix (not past the trimmed prefix)
    suffix = 0
    while (
        suffix < len(ref) - prefix
        and suffix < len(alt) - prefix
        and ref[len(ref) - 1 - suffix] == alt[len(alt) - 1 - suffix]
    ):
        suffix += 1
    ref_seg = ref[prefix : len(ref) - suffix]
    alt_seg = alt[prefix : len(alt) - suffix]
    return pos0 + prefix, ref_seg, alt_seg


def allele_depths(record, alt_index: int):
    """Return (ref_depth, alt_depth) from the first sample's AD field, or (None, None)"""
    if not len(record.samples):
        return None, None
    sample = next(iter(record.samples.values()))
    try:
        ad = sample["AD"]
    except (KeyError, ValueError, TypeError):
        ad = None
    if ad is None:
        return None, None
    ref_depth = ad[0] if len(ad) > 0 else None
    alt_depth = ad[alt_index + 1] if len(ad) > alt_index + 1 else None
    return ref_depth, alt_depth


def _depth_str(depth) -> str:
    return "NA" if depth is None else str(depth)


def genes_for_record(record, interval_index: dict) -> list:
    """Resolve the GeneModel(s) a record belongs to.
    Uses interval overlap"""
    models = []
    seen = set()
    for start, end, model in interval_index.get(record.contig, []):
        if record.start < end and record.stop > start and id(model) not in seen:
            seen.add(id(model))
            models.append(model)
    return models


def annotate_vcf(
    vcf, models_by_key: dict, suppress_downstream_frameshift: bool = True
) -> list:
    """Annotate every query-gene variant in a VCF.

    Returns a list of annotation dicts in VCF read order.  By default every
    variant that lies downstream (3') of a frameshift within the same gene is
    dropped, because a frameshift invalidates every codon after it; set
    ``suppress_downstream_frameshift=False`` to annotate those variants anyway."""
    interval_index = defaultdict(list)
    for model in set(models_by_key.values()):
        if model.genomic_start is not None:
            interval_index[model.contig].append(
                (model.genomic_start, model.genomic_end, model)
            )

    annotations = []
    read_order = 0
    for record in vcf:
        models = genes_for_record(record, interval_index)
        if not models or not record.alts or not is_nucleotide_allele(record.ref):
            continue
        for alt_index, alt in enumerate(record.alts):
            # skip symbolic ('<...>'), spanning-deletion ('*') and NON-REF alleles
            if not is_nucleotide_allele(alt):
                continue
            ref_depth, alt_depth = allele_depths(record, alt_index)
            changed_pos0, ref_seg, alt_seg = normalize_indel(
                record.start, record.ref, alt
            )
            for model in models:
                try:
                    ann = Variant(model, changed_pos0, ref_seg, alt_seg).annotate()
                except Exception as exc:  # never let one record abort the report
                    logger.warning(
                        f"Skipping {record.contig}:{record.pos} "
                        f"{record.ref}>{alt} for {model.gene_id}: {exc}"
                    )
                    ann = None
                if ann is None:
                    continue
                ann.update(
                    {
                        "gene_id": model.gene_id,
                        "product": model.product,
                        "ref_allele": record.ref,
                        "alt_allele": alt,
                        "ref_depth": ref_depth,
                        "alt_depth": alt_depth,
                        "read_order": read_order,
                    }
                )
                annotations.append(ann)
        read_order += 1

    if suppress_downstream_frameshift:
        return _apply_frameshift_suppression(annotations)
    annotations.sort(key=lambda a: a["read_order"])
    return annotations


def _apply_frameshift_suppression(annotations: list) -> list:
    """Drop variants that fall 3' of the first frameshift within each gene"""
    first_fs_pos = {}
    for ann in annotations:
        if ann["is_frameshift"]:
            gene = ann["gene_id"]
            pos = ann["cds_pos"]
            if gene not in first_fs_pos or pos < first_fs_pos[gene]:
                first_fs_pos[gene] = pos

    kept = []
    for ann in annotations:
        cutoff = first_fs_pos.get(ann["gene_id"])
        if cutoff is not None and ann["cds_pos"] > cutoff:
            logger.debug(
                f"Suppressing {ann['gene_id']} {ann['hgvs_c']} "
                f"(downstream of frameshift at c.{cutoff + 1})"
            )
            continue
        kept.append(ann)

    kept.sort(key=lambda a: a["read_order"])
    return kept


def format_report(annotations: list, ordered_query_genes) -> str:
    """Render the annotation report string (or the no-variants sentence)"""
    if not annotations:
        genes = ",".join(ordered_query_genes)
        return (
            f"No variants identified in queried genes ({genes}) "
            "relative to the reference genome"
        )

    entries = []
    for ann in annotations:
        depths = (
            f"{ann['ref_allele']}:{_depth_str(ann['ref_depth'])} "
            f"{ann['alt_allele']}:{_depth_str(ann['alt_depth'])}"
        )
        entries.append(
            f"{ann['gene_id']}: {ann['product']} "
            f"({ann['so_term']} {ann['hgvs_c']} {ann['hgvs_p']}; {depths})"
        )
    return ",".join(entries)


def run(
    vcffile: str,
    reference_gbff: str,
    reference_gff: str,
    reference_fa: str,
    query_genes,
    feature_qualifier: str = "product",
    exact_match: bool = False,
    transl_table: int = None,
    suppress_downstream_frameshift: bool = True,
    gene_vcf: str = "GENE_VARIANTS.vcf",
) -> str:
    """Annotate a VCF and return the report string.

    When ``gene_vcf`` is set (the default), the query-gene-overlapping variants
    are also written there with a ``GENE`` INFO field naming the overlapping
    gene(s)"""
    query_arg = query_genes if isinstance(query_genes, (list, tuple)) else [query_genes]
    ordered = ordered_query_genes(query_arg)
    # inefficiently reads VCF into memory
    vcf = pysam.VariantFile(vcffile)
    contig_names = set(vcf.header.contigs)
    if reference_gbff:
        models_by_key = build_gene_models_gbff(
            reference_gbff,
            contig_names,
            ordered,
            feature_qualifier,
            exact_match=exact_match,
            transl_table_override=transl_table,
        )
    elif reference_gff and reference_fa:
        models_by_key = build_gene_models_gff(
            reference_gff,
            reference_fa,
            contig_names,
            ordered,
            feature_qualifier,
            exact_match=exact_match,
            transl_table_override=transl_table,
        )
    else:
        raise FileNotFoundError("GBFF or GFF and FASTA not provided")
    # extract gene-overlapping variants to a VCF (offloaded from gene_coverage);
    # a model is registered under many keys, so dedupe to the distinct models
    if gene_vcf:
        extract_vcf_genes(
            vcffile, list(dict.fromkeys(models_by_key.values())), gene_vcf
        )
    annotations = annotate_vcf(
        vcf, models_by_key, suppress_downstream_frameshift=suppress_downstream_frameshift
    )
    return format_report(annotations, ordered)


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the variant_annotation arguments on ``parser``"""
    parser.add_argument("--vcf", required=True, help="VCF of gene-overlapping variants")
    parser.add_argument(
        "--reference_gbff",
        help="reference GenBank supplying coding sequence and coordinates",
    )
    parser.add_argument(
        "--reference_gff",
        help="reference GFF supplying coding coordinates"
    )
    parser.add_argument(
        "--reference_fa",
        help="reference genome FASTA"
    )
    parser.add_argument("--query_genes", nargs="+", required=True)
    parser.add_argument("--feature_qualifier", default="product")
    parser.add_argument("--exact_match", action="store_true")
    parser.add_argument(
        "--transl_table",
        type=int,
        default=None,
        help="override the genetic code (default: read /transl_table from each CDS)",
    )
    parser.add_argument(
        "--annotate_downstream_of_frameshift",
        action="store_true",
        help="annotate variants 3' of a frameshift (by default they are "
        "suppressed because a frameshift invalidates every downstream codon)",
    )
    parser.add_argument("--output", default="VARIANT_ANNOTATIONS.txt")
    parser.add_argument(
        "--gene_vcf",
        default="GENE_VARIANTS.vcf",
        help="write query-gene-overlapping variants (with a GENE INFO field) to "
        "this VCF (default: GENE_VARIANTS.vcf)",
    )
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Validate arguments, run the annotation and write/print the report"""
    has_gbff = bool(args.reference_gbff)
    has_gff_fa = bool(args.reference_gff and args.reference_fa)
    if not (has_gbff or has_gff_fa):
        raise ValueError("--reference_gbff OR --reference_gff and --reference_fa required")
    elif has_gbff and (args.reference_gff or args.reference_fa):
        raise ValueError("--reference_gbff is mutually exclusive with --reference_gff and --reference_fa")

    report = run(
        args.vcf,
        args.reference_gbff,
        args.reference_gff,
        args.reference_fa,
        args.query_genes,
        feature_qualifier=args.feature_qualifier,
        exact_match=args.exact_match,
        transl_table=args.transl_table,
        suppress_downstream_frameshift=not args.annotate_downstream_of_frameshift,
        gene_vcf=args.gene_vcf,
    )
    with open(args.output, "w") as out:
        out.write(report + "\n")
    print(report)
    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.variant_annotation``)"""
    parser = argparse.ArgumentParser(
        description="Annotate protein-level consequences of gene-overlapping variants "
        "(task_variant_annotate.wdl dependency; github.com/theiagen/public_health_bioinformatics)"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
