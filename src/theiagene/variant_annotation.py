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
from theiagene.lib.vcf import extract_vcf_genes, flatten_coords_by_contig  # noqa: F401
from theiagene.lib.gff import iter_gff_features
from theiagene.lib.logging_config import configure_logging

# domain classes, now living in lib; re-exported so callers/tests reach them here
from theiagene.lib.gene_model import GeneModel  # noqa: F401
from theiagene.lib.variant import Variant


logger = logging.getLogger(__name__)


def _ordered_genomic_positions(parts, strand):
    """Ordered 0-based genomic positions of a CDS in translation (5'->3') order"""
    ordered = sorted(((int(s), int(e)) for s, e in parts), key=lambda x: x[0])
    positions = []
    for start, end in ordered:
        positions.extend(range(start, end))
    if strand == -1:
        positions.reverse()
    return positions


def build_gene_models_gbff(
    reference_gbff: str,
    contig_names: set,
    query_genes,
    feature_type: str,
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
    id_qualifiers = [feature_qualifier.strip(), "gene", "locus_tag", "protein_id"]
    models_by_key = {}
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_seq = str(record.seq)
            for feature in record.features:
                if feature.type.lower() != feature_type.lower():
                    continue
                identifiers = []
                for qualifier in id_qualifiers:
                    identifiers.extend(feature.qualifiers.get(qualifier, []))
                if not identifiers:
                    continue

                matched_query = match_query(query_list, identifiers, exact_match)
                if matched_query is None:
                    continue

                strand = feature.location.strand
                if strand not in (1, -1):
                    logger.warning(
                        f"Skipping '{matched_query}' on {record.name}: "
                        f"unresolved strand ({strand})"
                    )
                    continue

                if transl_table_override is not None:
                    transl_table = transl_table_override
                else:
                    transl_table = int(
                        feature.qualifiers.get("transl_table", ["1"])[0]
                    )

                product_vals = feature.qualifiers.get(feature_qualifier.strip())
                product = product_vals[0] if product_vals else matched_query
                gene_id = normalize_name(matched_query)

                record_id = record.id
                # check appropriate query is use for VCF contigs
                if record_id not in contig_names:
                    record_id = record.name
                    if record_id not in contig_names:
                        raise KeyError(f"{record.id} and {record.name} not in VCF")

                model = GeneModel(
                    gene_id=gene_id,
                    product=product,
                    contig=record_id,
                    strand=strand,
                    transl_table=transl_table,
                )
                parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
                model.genomic_positions = _ordered_genomic_positions(parts, strand)
                model.finalize(contig_seq)

                keys = set()
                for ident in identifiers + [matched_query, product, gene_id]:
                    keys.update(
                        (ident, sanitize_info_value(ident), normalize_name(ident))
                    )
                for key in keys:
                    if not key:
                        continue
                    if key in models_by_key and models_by_key[key] is not model:
                        logger.warning(
                            f"'{key}' recovered multiple times; keeping first"
                        )
                    else:
                        models_by_key[key] = model
    return models_by_key


# GFF3 attribute keys used to group multi-segment CDS lines into one gene, in
# preference order; all segments of a CDS share these (Parent/ID especially)
_GFF_GROUP_KEYS = ("Parent", "ID", "locus_tag", "gene", "protein_id")


def _group_gff_cds(reference_gff: str, feature_type: str, id_qualifiers) -> list:
    """Group GFF CDS lines into one entry per gene.

    A multi-exon CDS is spread over several GFF lines that share an ``ID`` (or
    ``Parent``); this collapses them into a single group carrying every segment
    and the identifiers used for query matching.  Returns the groups in
    first-seen order for deterministic model registration."""
    groups = {}
    order = []
    for feature in iter_gff_features(reference_gff, feature_type):
        attrs = feature["attributes"]
        identifiers = [attrs[q] for q in id_qualifiers if q in attrs]
        if not identifiers:
            continue
        # pick the first available grouping key so exon segments coalesce
        group_id = None
        for key in _GFF_GROUP_KEYS:
            if key in attrs:
                group_id = (feature["seqid"], key, attrs[key])
                break
        if group_id is None:
            # no grouping attribute: treat this single line as its own gene
            group_id = (feature["seqid"], "line", feature["start"], feature["end"])
        group = groups.get(group_id)
        if group is None:
            group = {
                "seqid": feature["seqid"],
                "strand": feature["strand"],
                "attributes": attrs,
                "identifiers": [],
                "parts": [],
            }
            groups[group_id] = group
            order.append(group_id)
        group["parts"].append((feature["start"], feature["end"]))
        for ident in identifiers:
            if ident not in group["identifiers"]:
                group["identifiers"].append(ident)
    return [groups[group_id] for group_id in order]


def build_gene_models_gff(
    reference_gff: str,
    reference_fa: str,
    contig_names: set,
    query_genes,
    feature_type: str,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GFF and genome FA into {lookup_key: GeneModel} for every query gene.

    A CDS is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` attributes, so query sets may mix product
    names and locus tags.  Multi-exon CDS lines are grouped (by ``ID``/``Parent``)
    and assembled in translation order into a single coding model, mirroring the
    GBFF backend.  A model is registered under many lookup keys (raw, sanitized
    and normalized forms of every identifier) so it can be recovered from
    whatever identifier the VCF ``GENE`` field carries.

    The coding phase column is not applied; as with the GBFF backend, each CDS
    is assumed to begin on a codon boundary."""
    query_list = list(query_genes)
    qualifier = feature_qualifier.strip()
    # this may need to be exposed to enable user-modification
    id_qualifiers = [qualifier, "gene", "locus_tag", "protein_id"]
    fa_dict = SeqIO.to_dict(SeqIO.parse(reference_fa, "fasta"))
    models_by_key = {}

    for group in _group_gff_cds(reference_gff, feature_type, id_qualifiers):
        identifiers = group["identifiers"]
        matched_query = match_query(query_list, identifiers, exact_match)
        if matched_query is None:
            continue

        record_id = group["seqid"]
        strand = group["strand"]
        if strand not in (1, -1):
            logger.warning(
                f"Skipping '{matched_query}' on {record_id}: "
                f"unresolved strand ({strand})"
            )
            continue

        # check appropriate contig is used for the VCF and available in the FASTA
        if record_id not in contig_names:
            raise KeyError(f"{record_id} not in VCF")
        if record_id not in fa_dict:
            raise KeyError(f"{record_id} not in reference FASTA")
        contig_seq = str(fa_dict[record_id].seq)

        attrs = group["attributes"]
        if transl_table_override is not None:
            transl_table = transl_table_override
        else:
            transl_table = int(attrs.get("transl_table", "1"))

        product = attrs.get(qualifier, matched_query)
        gene_id = normalize_name(matched_query)

        model = GeneModel(
            gene_id=gene_id,
            product=product,
            contig=record_id,
            strand=strand,
            transl_table=transl_table,
        )
        model.genomic_positions = _ordered_genomic_positions(group["parts"], strand)
        model.finalize(contig_seq)

        keys = set()
        for ident in identifiers + [matched_query, product, gene_id]:
            keys.update(
                (ident, sanitize_info_value(ident), normalize_name(ident))
            )
        for key in keys:
            if not key:
                continue
            if key in models_by_key and models_by_key[key] is not model:
                logger.warning(
                    f"'{key}' recovered multiple times; keeping first"
                )
            else:
                models_by_key[key] = model
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
    feature_type: str = "CDS",
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
            feature_type,
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
            feature_type,
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
    parser.add_argument("--feature_type", default="CDS")
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
        feature_type=args.feature_type,
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
