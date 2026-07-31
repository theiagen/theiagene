"""Render VEP variant annotations into product-named report lines.

Given a VEP ``--tab`` output TSV and the reference GFF used to produce it, this
command keeps the rows whose ``Consequence`` is not suppressed, resolves each
row's ``Feature`` (an RNA) back to its CDS product name through the GFF feature
hierarchy, and rewrites the transcript/protein prefixes of the HGVSc/HGVSp
strings to that product name, e.g.::

    lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg)

When the source VCF is supplied via ``--vcf``, each line also carries
the variant's per-allele read depths, e.g.::

    lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)

Rows carrying neither an HGVSc nor an HGVSp string are ignored for now."""

import sys
import logging
import argparse

from theiagene.lib.feature import FeatureCol
from theiagene.lib.parsers import assimilate_gff, import_vcf
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)

# VEP's undefined-column placeholder
_UNDEFINED = {"-", "", "."}


def _split_qualifiers(raw: str) -> list:
    """Split a comma-/space-delimited qualifier string into individual keys,
    dropping empty tokens."""
    return raw.replace(",", " ").split()


def parse_vep_tsv(vep_tsv: str):
    """Yield each VEP ``--tab`` data row as a ``{column: value}`` dict.

    VEP emits ``##`` meta lines, then a single ``#``-prefixed column-header line,
    then tab-delimited data rows; the header names the columns the data rows are
    zipped against."""
    columns = None
    with open(vep_tsv) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                # the lone '#'-prefixed line carries the column names
                columns = line[1:].split("\t")
                continue
            if columns is None:
                raise ValueError(f"no column header found before data in {vep_tsv}")
            yield dict(zip(columns, line.split("\t")))


def _allele_depths(sample):
    """Return ``(ref_depth, [alt_depth, ...])`` for a VCF record's single sample,
    or None when the record carries no per-allele read counts.

    FreeBayes writes these as the ``AD`` FORMAT field (``ref, alt1, alt2, ...``);
    if ``AD`` is absent it falls back to the ``RO`` (reference-observation) and
    ``AO`` (alternate-observation) fields. gVCF reference blocks carry neither and
    yield None."""
    ad = sample.get("AD")
    if ad:
        return ad[0], list(ad[1:])
    ro = sample.get("RO")
    ao = sample.get("AO")
    if ro is not None and ao is not None:
        ao = ao if isinstance(ao, (tuple, list)) else (ao,)
        return ro, list(ao)
    return None


def build_depth_index(vcf: str) -> dict:
    """Map each VCF record's raw ``(contig, pos, ref)`` to its per-allele depths.

    The value is ``(ref_depth, {alt: alt_depth})``. Keying on the *raw*
    (unnormalized) coordinates is deliberate: VEP echoes exactly these back in its
    ``Uploaded_variation`` column (``contig_pos_ref/alt``), so a report row can be
    matched to its source record without redoing VEP's allele normalization.
    Records without per-allele depths (e.g. gVCF reference blocks) are skipped."""
    index = {}
    handle = import_vcf(vcf)
    try:
        for record in handle:
            sample = next(iter(record.samples.values()), None)
            if sample is None:
                continue
            depths = _allele_depths(sample)
            if depths is None:
                continue
            ref_depth, alt_depths = depths
            alts = record.alts or ()
            index[(record.contig, record.pos, record.ref)] = (
                ref_depth,
                dict(zip(alts, alt_depths)),
            )
    finally:
        handle.close()
    return index


def _depth_suffix(row: dict, depth_index: dict):
    """Render ``"{ref}:{ref_depth} {alt}:{alt_depth}"`` for a VEP row, or None when
    the row's variant is absent from ``depth_index`` or its allele can't be
    resolved.

    ``Uploaded_variation`` is ``contig_pos_ref/alt[,alt...]`` in VEP's raw input
    coordinates; the contig/pos/ref triple keys the index and the row's ``Allele``
    (falling back to the sole alt of a biallelic site) picks which alt's depth to
    report."""
    uploaded = row.get("Uploaded_variation", "")
    try:
        rest, _alts = uploaded.rsplit("/", 1)
        prefix, ref = rest.rsplit("_", 1)
        contig, pos = prefix.rsplit("_", 1)
        pos = int(pos)
    except ValueError:
        return None
    entry = depth_index.get((contig, pos, ref))
    if entry is None:
        return None
    ref_depth, alt_depths = entry
    allele = row.get("Allele")
    if allele in alt_depths:
        alt = allele
    elif len(alt_depths) == 1:
        alt = next(iter(alt_depths))
    else:
        # can't disambiguate which alt of a multiallelic site this row is about
        return None
    return f"{ref}:{ref_depth} {alt}:{alt_depths[alt]}"


def _consequences(row: dict) -> list:
    """Split a row's ``Consequence`` column into its individual terms (VEP joins
    co-occurring consequences with ',' or '&')."""
    raw = row.get("Consequence", "")
    return raw.replace("&", ",").replace(" ", "").split(",") if raw else []


def _cds_product(features: FeatureCol, feature_id: str, feature_type: str, qualifiers: list):
    """Resolve a VEP ``Feature`` (an RNA id) to the ``feature_qualifier`` value on
    one of its ``feature_type`` (CDS) descendants.

    The attribute key is matched case-insensitively; the first matching value on
    the first qualifying descendant wins. Returns None if the feature is absent
    from the collection or carries no such descendant/qualifier."""
    feature = features.get(feature_id)
    if feature is None:
        return None
    wanted = {qualifier.lower() for qualifier in qualifiers}
    # only the requested-type subfeatures beneath this RNA (e.g. CDS)
    for subfeature in FeatureCol(feature.descendants, group=False)[feature_type]:
        for key, value in subfeature.attributes.items():
            if value and key.lower() in wanted:
                return value
    return None


def _hgvs_suffix(value: str, strip_parens: bool = False):
    """Return the portion of an HGVS string after the ``transcript:`` /
    ``protein:`` prefix, or None when the column is undefined.

    ``strip_parens`` drops the parentheses VEP wraps predicted protein changes in
    (``p.(Lys143Arg)`` -> ``p.Lys143Arg``)."""
    if value is None or value in _UNDEFINED:
        return None
    suffix = value.split(":", 1)[1] if ":" in value else value
    if strip_parens:
        suffix = suffix.replace("(", "").replace(")", "")
    return suffix


def report_line(row: dict, product: str, depth_index: dict = None) -> str:
    """Build the product-named report line for a single kept VEP row.

    The transcript/protein prefixes of HGVSc/HGVSp are translated to ``product``
    (spaces -> '.') and factored out as a shared label; the consequence and the
    surviving HGVS suffixes follow in parentheses. When ``depth_index`` is given
    and the row's variant resolves in it, the per-allele read depths are appended
    after a ``;`` (e.g. ``... p.Lys143Arg; T:0 C:562``)."""
    dotted = product.replace(" ", ".")
    consequence = row.get("Consequence", "")
    pieces = [consequence]
    hgvsc = _hgvs_suffix(row.get("HGVSc"))
    hgvsp = _hgvs_suffix(row.get("HGVSp"), strip_parens=True)
    if hgvsc:
        pieces.append(hgvsc)
    if hgvsp:
        pieces.append(hgvsp)
    body = " ".join(pieces)
    if depth_index is not None:
        depths = _depth_suffix(row, depth_index)
        if depths:
            body += f"; {depths}"
    return f"{dotted}: {product} ({body})"


def report_variants(
    vep_tsv: str,
    features: FeatureCol,
    suppress: set,
    feature_type: str,
    qualifiers: list,
    depth_index: dict = None,
) -> list:
    """Turn a VEP TSV into product-named report lines.

    A row is dropped when any of its consequence terms is suppressed, when it
    carries neither an HGVSc nor an HGVSp string, or when its ``Feature`` cannot
    be resolved to a CDS product in ``features``. When ``depth_index`` is given
    (see :func:`build_depth_index`), each kept line carries the variant's
    per-allele read depths."""
    lines = []
    for row in parse_vep_tsv(vep_tsv):
        if any(consequence in suppress for consequence in _consequences(row)):
            continue
        # nothing to translate without at least one HGVS string
        if row.get("HGVSc") in _UNDEFINED and row.get("HGVSp") in _UNDEFINED:
            continue
        feature_id = row.get("Feature")
        product = _cds_product(features, feature_id, feature_type, qualifiers)
        if product is None:
            logger.warning(
                f"no {feature_type} {qualifiers} qualifier resolved for feature "
                f"{feature_id!r}; skipping variant {row.get('Uploaded_variation')}"
            )
            continue
        lines.append(report_line(row, product, depth_index))
    return lines


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the report_variants arguments on ``parser``"""
    parser.add_argument("--vep_tsv", required=True)
    parser.add_argument("--reference_gff", required=True)
    parser.add_argument(
        "--vcf",
        help="VCF the VEP annotations were called from; when given, each report "
        "line carries the variant's per-allele read depths (e.g. 'T:0 C:562')",
    )
    parser.add_argument(
        "--suppress",
        default="intergenic_variant,coding_sequence_variant",
        help="comma-/space-delimited consequence type(s) whose rows are dropped",
    )
    parser.add_argument(
        "--feature_qualifier",
        default="product",
        help="attribute key(s), matched case-insensitively, read off the CDS "
        "descendant(s) to name the product",
    )
    parser.add_argument("--feature_type", default="CDS")
    parser.add_argument(
        "--output",
        help="write report lines here instead of stdout",
    )
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Render the VEP annotations into product-named report lines"""
    features = assimilate_gff(args.reference_gff)
    suppress = set(_split_qualifiers(args.suppress))
    qualifiers = _split_qualifiers(args.feature_qualifier)
    depth_index = build_depth_index(args.vcf) if args.vcf else None

    lines = report_variants(
        args.vep_tsv, features, suppress, args.feature_type, qualifiers, depth_index
    )

    if args.output:
        with open(args.output, "w") as handle:
            handle.write("\n".join(lines) + ("\n" if lines else ""))
        logger.debug(f"Wrote {len(lines)} variant report line(s) to {args.output}")
    else:
        for line in lines:
            print(line)
        logger.debug(f"Reported {len(lines)} variant(s)")

    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.report_variants``)"""
    parser = argparse.ArgumentParser(
        description="Extracts variant annotations from VEP output report"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
