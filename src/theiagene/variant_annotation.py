"""Annotate the protein-level consequences of variants that overlap query genes.

This is the ``variant_annotation`` subcommand of the ``theiagene`` entrypoint
(formerly the standalone ``variant_annotation.py`` script).  Given a VCF, a
reference GenBank (GBFF) that supplies the coding sequence, strand, product name
and translation table for each gene, and a set of query genes, it reports the
effect of every variant on the translated protein.

Query genes are resolved for each variant from the VCF ``GENE`` INFO field when
one is present (for example, a VCF pre-filtered by ``gene_coverage``) and
otherwise by interval overlap against the coding models built from the GBFF, so
it works equally on a pre-extracted or a raw VCF.

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
from Bio import SeqIO

# shared helpers, re-exported so callers (and tests) can reach them here
from theiagene.lib.sequence import (  # noqa: F401
    is_nucleotide_allele,
    aa3,
    complement,
    translate,
)
from theiagene.lib.query import (  # noqa: F401
    normalize_name,
    sanitize_info_value,
    match_query,
    ordered_query_genes,
)
from theiagene.lib.parsers import (  # noqa: F401
    iter_gbff_raw,
    iter_gff_raw,
    match_identifiers,
    resolve_contig,
    extract_vcf_genes,
    flatten_coords_by_contig,
)
from theiagene.lib.logging_config import configure_logging

# domain classes, now living in lib; re-exported so callers/tests reach them here
from theiagene.lib.gene_model import GeneModel  # noqa: F401
from theiagene.lib.variant import Variant


logger = logging.getLogger(__name__)


def _id_qualifiers(feature_qualifier: str) -> list:
    """CDS qualifiers matched for query resolution (query sets may mix product
    names and locus tags), the feature qualifier first"""
    return [feature_qualifier.strip(), "gene", "locus_tag", "protein_id"]


def _register_model(models_by_key: dict, model: GeneModel, key_sources) -> None:
    """Register a model under raw/sanitized/normalized forms of every identifier,
    so it can be recovered from whatever identifier the VCF ``GENE`` field carries"""
    keys = set()
    for ident in key_sources:
        keys.update((ident, sanitize_info_value(ident), normalize_name(ident)))
    for key in keys:
        if not key:
            continue
        if key in models_by_key and models_by_key[key] is not model:
            logger.warning(f"'{key}' recovered multiple times; keeping first")
        else:
            models_by_key[key] = model


def _assemble_model(raw, contig, matched_query, feature_qualifier, transl_table_override):
    """Resolve a matched Gene's identity and generate its GeneModel.
    Returns (model, gene_id, product)"""
    qualifier = feature_qualifier.strip()
    if transl_table_override is not None:
        transl_table = transl_table_override
    else:
        transl_table = int(raw.qualifiers.get("transl_table", ["1"])[0])
    product_vals = raw.qualifiers.get(qualifier)
    product = product_vals[0] if product_vals else matched_query
    gene_id = normalize_name(matched_query)
    # stamp the resolved identity onto the parsed Gene, then let the GeneModel
    # derive itself from it (it carries the CDS coordinates and reference sequence)
    raw.gene_id = gene_id
    raw.product = product
    raw.contig = contig
    raw.transl_table = transl_table
    model = GeneModel.from_gene(raw)
    return model, gene_id, product


def build_gene_models_gbff(
    reference_gbff: str,
    contig_names: set,
    query_genes,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GBFF into {lookup_key: GeneModel} for every query gene.

    A CDS is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` qualifiers, so query sets may mix product
    names and locus tags.  A model is registered under many lookup keys (raw,
    sanitized and normalized forms of every identifier) so it can be recovered
    from whatever identifier the VCF ``GENE`` field carries."""
    query_list = list(query_genes)
    id_qualifiers = _id_qualifiers(feature_qualifier)
    models_by_key = {}
    for raw in iter_gbff_raw(reference_gbff):
        matched_query, identifiers = match_identifiers(
            raw.qualifiers, query_list, id_qualifiers,
            exact_match=exact_match, normalize=True,
        )
        if matched_query is None:
            continue
        if raw.strand not in (1, -1):
            logger.warning(
                f"Skipping '{matched_query}': unresolved strand ({raw.strand})"
            )
            continue
        if not raw.cds:
            logger.warning(f"Skipping '{matched_query}': no CDS coordinates to model")
            continue
        contig = resolve_contig(raw.contig_candidates, contig_names, True, "VCF")
        model, gene_id, product = _assemble_model(
            raw, contig, matched_query, feature_qualifier, transl_table_override
        )
        _register_model(
            models_by_key, model, identifiers + [matched_query, product, gene_id]
        )
    return models_by_key


def build_gene_models_gff(
    reference_gff: str,
    reference_fa: str,
    contig_names: set,
    query_genes,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GFF and genome FA into {lookup_key: GeneModel} for every query gene.

    A CDS is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` attributes, so query sets may mix product
    names and locus tags.  Multi-exon CDS are assimilated through the GFF3
    parent/child hierarchy (``CDS -> RNA -> gene``) and assembled in translation
    order into a single coding model, mirroring the GBFF backend.  A model is
    registered under many lookup keys (raw, sanitized and normalized forms of
    every identifier) so it can be recovered from whatever identifier the VCF
    ``GENE`` field carries.

    The coding phase column is not applied; as with the GBFF backend, each CDS
    is assumed to begin on a codon boundary."""
    query_list = list(query_genes)
    id_qualifiers = _id_qualifiers(feature_qualifier)
    fa_dict = SeqIO.to_dict(SeqIO.parse(reference_fa, "fasta"))
    models_by_key = {}
    for raw in iter_gff_raw(reference_gff, fa_dict):
        matched_query, identifiers = match_identifiers(
            raw.qualifiers, query_list, id_qualifiers,
            exact_match=exact_match, normalize=True,
        )
        if matched_query is None:
            continue
        if raw.strand not in (1, -1):
            logger.warning(
                f"Skipping '{matched_query}': unresolved strand ({raw.strand})"
            )
            continue
        if not raw.cds:
            logger.warning(f"Skipping '{matched_query}': no CDS coordinates to model")
            continue
        # check appropriate contig is used for the VCF and available in the FASTA
        contig = resolve_contig(raw.contig_candidates, contig_names, True, "VCF")
        if raw.contig_seq is None:
            raise KeyError(f"{contig} not in reference FASTA")
        model, gene_id, product = _assemble_model(
            raw, contig, matched_query, feature_qualifier, transl_table_override
        )
        _register_model(
            models_by_key, model, identifiers + [matched_query, product, gene_id]
        )
    return models_by_key


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
) -> str:
    """Annotate a VCF and return the report string"""
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
    )
    with open(args.output, "w") as out:
        out.write(report + "\n")
    print(report)
    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.variant_annotation``)"""
    parser = argparse.ArgumentParser(
        description="Annotate protein-level consequences of gene-overlapping variants "
        "(task_gene_coverage.wdl dependency; github.com/theiagen/public_health_bioinformatics)"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
