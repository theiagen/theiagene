"""Render VEP variant annotations into gene-labelled report lines.

Given a VEP ``--tab`` output TSV and the reference GFF used to produce it, this
command keeps the rows whose ``Consequence`` is not suppressed, resolves each
row's ``Location`` back to a CDS product name through the GFF feature hierarchy,
and replaces the transcript/protein prefixes of the HGVSc/HGVSp strings with a
single leading gene label, e.g.::

    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg)

The label is the ``--query_genes`` term that matches the row's feature -- the
name the caller asked about, not the full product it resolved to -- and the
product follows in quotes so a name carrying commas survives being joined into
a comma-delimited report. Without ``--query_genes`` (or when none match) the
label falls back to the normalized product name.

A row is tied back to the annotation by its ``Location`` rather than by its
``Feature`` column: the identifier VEP writes there is whichever attribute its
own GFF parser read off the transcript, which varies by annotation source and
need not be the ``ID`` this package keys features on. ``Feature`` is consulted
only to break a tie between several annotation units overlapping one variant.

When the source VCF is supplied via ``--vcf``, each line also carries
the variant's per-allele read depths, e.g.::

    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)

Rows carrying neither an HGVSc nor an HGVSp string are ignored for now."""

import re
import sys
import logging
import argparse
from urllib.parse import unquote

from theiagene.lib.feature import FeatureCol
from theiagene.lib.parsers import assimilate_gff, import_vcf
from theiagene.lib.query import (
    ordered_query_genes,
    extract_queries_from_bed,
    feature_identifiers,
    iter_descendants,
    match_query,
    normalize_name,
    split_qualifiers,
)
from theiagene.lib.logging_config import configure_logging


logger = logging.getLogger(__name__)

# VEP's undefined-column placeholders. None is included so a column that is
# absent from the TSV entirely (VEP run without --hgvs) or missing from a short
# data row -- both of which surface as a `.get` miss -- reads as undefined
# rather than as a value, which would otherwise let a row carrying no HGVS
# string at all past the filter in `report_variants`.
_UNDEFINED = {"-", "", ".", None}

# the attribute keys VEP's `Feature` column can be drawn from; which one its GFF
# parser lands on varies by annotation source, so each is tried in turn
_UNIT_ID_KEYS = ("ID", "Name", "transcript_id")


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


def _type_qualifier(feature, feature_type: str, qualifiers: list):
    """Resolve an annotation unit to the ``feature_qualifier`` value carried by
    one of its ``feature_type`` (CDS) records.

    The unit is whichever level the annotation resolved the variant to -- an RNA,
    a gene, or a bare ``feature_type`` record -- so the search covers the unit
    itself and every descendant beneath it at any depth: the unit itself so a
    bare CDS carrying no gene parent resolves its own product, and the recursive
    walk so a gene sitting over gene -> RNA -> CDS still reaches the CDS. This is
    the same set of subfeatures :func:`theiagene.lib.query._grouped_query_ranges`
    reads coordinates from, so the report and the extraction resolve a unit
    through one walk.

    The attribute key is matched case-insensitively; the first matching value on
    the first qualifying record wins. Returns None if the feature is None or
    carries no such record/qualifier."""
    if feature is None:
        return None
    wanted = {qualifier.lower() for qualifier in qualifiers}
    # the requested-type records at or beneath this unit (e.g. CDS); the class
    # filter drops the unit itself unless it is already of that type
    subfeatures = FeatureCol([feature] + list(iter_descendants(feature)), group=False)
    for subfeature in subfeatures[feature_type]:
        for key, value in subfeature.attributes.items():
            if value and key.lower() in wanted:
                return value
    return None


def _query_label(feature, query_list: list, qualifiers: list, exact_match: bool):
    """Return the ``--query_genes`` term that selects ``feature``, or None.

    The term is matched against every identifier the feature, its parent and its
    descendants carry (the same candidates :mod:`theiagene.lib.query` matches
    when extracting the variants), so the report is labelled by the name the
    user asked for -- 'FKS1' rather than the full CDS product it resolved to.
    None is returned when no queries were supplied or none of them match, which
    leaves the caller to fall back to the product-derived label."""
    if not query_list or feature is None:
        return None
    return match_query(query_list, feature_identifiers(feature, qualifiers), exact_match)


# HGVS writes a synonymous protein change as ``p.Asp164=`` ('=' meaning "residue
# unchanged"); this captures the prefix, the reference residue (3- or 1-letter
# code) and its position so the residue can be repeated in place of the '='
_SYNONYMOUS = re.compile(r"^(p\.)([A-Z][a-z]{2}|[A-Z])(\d+)=$")


def _expand_synonymous(suffix: str) -> str:
    """Rewrite a synonymous HGVS protein change into its expanded form
    (``p.Asp164=`` -> ``p.Asp164Asp``).

    A change naming no single reference residue -- ``p.=`` (whole protein
    unchanged) or a range such as ``p.Asp164_Leu166=`` -- has nothing
    unambiguous to repeat, so it is returned untouched. Anything that is not a
    synonymous change passes through unchanged."""
    match = _SYNONYMOUS.match(suffix)
    if not match:
        return suffix
    prefix, residue, position = match.groups()
    return f"{prefix}{residue}{position}{residue}"


def _hgvs_suffix(value: str, strip_parens: bool = False):
    """Return the portion of an HGVS string after the ``transcript:`` /
    ``protein:`` prefix, or None when the column is undefined.

    VEP percent-encodes the characters that are reserved in a VCF INFO field
    (``=``, ``;``, ``,``, ``&``, ``%``) even in its ``--tab`` output, so the
    suffix is URL-decoded here -- otherwise a synonymous change arrives as
    ``p.Asp164%3D``. The decoded ``=`` is then expanded to repeat the reference
    residue (see :func:`_expand_synonymous`).

    ``strip_parens`` drops the parentheses VEP wraps predicted protein changes in
    (``p.(Lys143Arg)`` -> ``p.Lys143Arg``)."""
    if value in _UNDEFINED:
        return None
    suffix = unquote(value.split(":", 1)[1] if ":" in value else value)
    if strip_parens:
        suffix = suffix.replace("(", "").replace(")", "")
    return _expand_synonymous(suffix)


def report_line(row: dict, product: str, depth_index: dict = None, label: str = None) -> str:
    """Build the report line for a single kept VEP row.

    The line leads with ``label`` -- the query-gene term the variant was
    extracted under (e.g. ``FKS1``) -- followed by the quoted CDS product; the
    consequence and the surviving HGVS suffixes follow in parentheses. Quoting
    the product keeps a name carrying commas ('1,3-beta-glucan synthase ...')
    intact once the caller joins the lines with ','. Without a ``label`` (no
    ``--query_genes``, or none matched) the product itself is normalized into
    one (spaces -> '.'), preserving the previous behaviour. When ``depth_index``
    is given and the row's variant resolves in it, the per-allele read depths
    are appended after a ``;`` (e.g. ``... p.Lys143Arg; T:0 C:562``)."""
    label = label or normalize_name(product)
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
    return f'{label}: "{product}" ({body})'


def _extract_location(row: dict):
    """Parse VEP's ``Location`` into ``(contig, start, end)``, 0-based half-open.

    VEP writes ``contig:pos`` for a substitution, ``contig:start-end`` for a
    deletion, and ``contig:start-end`` with start > end for an insertion (the
    zero-length span between two bases). Coordinates are 1-based and inclusive,
    so only the start shifts; sorting the pair turns an insertion into the
    two-base window it sits between, which is what an overlap test needs. The
    contig is split from the right so a name containing ':' survives.

    None is returned for a ``Location`` that is undefined, absent from the TSV
    entirely, or unparseable, leaving the caller to skip the row rather than
    resolve it against a bogus coordinate."""
    raw = row.get("Location")
    if raw in _UNDEFINED:
        return None
    contig, sep, span = raw.rpartition(":")
    if not sep or not contig:
        return None
    low, _, high = span.partition("-")
    try:
        low = int(low)
        # a bare position names a single base, so both ends come from it
        high = int(high) if high else low
    except ValueError:
        return None
    return contig, min(low, high) - 1, max(low, high)


def _select_unit(units: list, vep_feature: str):
    """Pick the single annotation unit a VEP row is about, or None.

    A row's HGVS strings are computed against exactly one transcript, and VEP has
    already emitted a separate row per transcript, so a row must resolve to one
    unit rather than fan out over every unit the variant overlaps. A sole
    overlapping unit is taken as-is -- the coordinate has already done the work
    and there is nothing to disambiguate. Only where several overlap does the
    row's ``Feature`` break the tie, compared against the unit's own identifiers
    (never its parent's or its CDS product's, which sibling transcripts of one
    gene share and which would therefore match both).

    None is returned when nothing overlapped, when an ambiguous overlap carries
    no ``Feature`` to arbitrate it, and when none of the units answer to the one
    it carries -- all of which leave the caller to skip the row rather than
    attribute it to whichever unit happened to sort first."""
    if len(units) < 2:
        return units[0] if units else None
    if vep_feature in _UNDEFINED:
        return None
    for unit in units:
        # `fid` alongside the attributes: group_features may have renamed a
        # colliding ID, leaving the one VEP saw only in `attributes`
        identifiers = [unit.fid] + [unit.attributes.get(key) for key in _UNIT_ID_KEYS]
        if any(identifier == vep_feature for identifier in identifiers):
            return unit
    return None


def report_variants(
    vep_tsv: str,
    features: FeatureCol,
    suppress: set,
    feature_type: str,
    qualifiers: list,
    depth_index: dict = None,
    query_list: list = None,
    exact_match: bool = False,
) -> list:
    """Turn a VEP TSV into query-labelled report lines.

    A row is dropped when any of its consequence terms is suppressed, when it
    carries neither an HGVSc nor an HGVSp string, when its ``Location`` does not
    parse, or when that location cannot be resolved to a CDS product in
    ``features`` -- either because no annotation unit overlaps it, because
    several do and none answers to the row's ``Feature`` (see
    :func:`_select_unit`), or because the resolved unit carries no qualifier.
    Each kept line leads with the ``query_list`` term that matched the row's
    unit, falling back to the product-derived label when no query matched. When
    ``depth_index`` is given (see :func:`build_depth_index`), each kept line
    carries the variant's per-allele read depths."""
    lines = []
    for row in parse_vep_tsv(vep_tsv):
        if any(consequence in suppress for consequence in _consequences(row)):
            continue
        # nothing to translate without at least one HGVS string
        if row.get("HGVSc") in _UNDEFINED and row.get("HGVSp") in _UNDEFINED:
            continue
        # the variant is placed by coordinate rather than by the row's `Feature`,
        # whose spelling depends on which attribute VEP's GFF parser read
        location = _extract_location(row)
        if location is None:
            logger.warning(f"cannot resolve location {row.get('Location')}; "
                           f"skipping variant {row.get('Uploaded_variation')}")
            continue
        seqid, start, end = location
        hits = features.index(seqid, start, end)
        # `index` returns every overlapping record regardless of class, so units
        # are taken from one level: the transcripts VEP itself annotates against,
        # else the genes of an annotation carrying no transcript level (gene ->
        # CDS), else bare feature_type records owning no gene at all. Taking one
        # level keeps a locus from counting once per tier of its own hierarchy
        units = hits.rnas or hits.genes or hits[feature_type]
        feature = _select_unit(units, row.get("Feature"))
        qualifier_hit = _type_qualifier(feature, feature_type, qualifiers)
        if qualifier_hit is None:
            logger.warning(
                f"no {feature_type} {qualifiers} qualifier resolved for feature "
                f"{row.get('Location')}; skipping variant {row.get('Uploaded_variation')}"
            )
            continue
        label = _query_label(feature, query_list, qualifiers, exact_match)
        if query_list and label is None:
            logger.warning(
                f"no query gene matched feature {row.get('Location')}; labelling variant "
                f"{row.get('Uploaded_variation')} by its product instead"
            )
        lines.append(report_line(row, qualifier_hit, depth_index, label))
    return lines


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the report_variants arguments on ``parser``"""
    parser.add_argument("--vep_tsv", required=True)
    parser.add_argument("--reference_gff", required=True)
    parser.add_argument(
        "--query_genes",
        help="comma-delimited query gene name(s)
    )
    parser.add_argument(
        "--bedfile",
        help="BED whose name column supplies the query names when --query_genes "
        "is not given",
    )
    parser.add_argument(
        "--exact_match",
        action="store_true",
        help="require a query term to equal a feature identifier rather than be "
        "a substring of one (always case-sensitive)",
    )
    parser.add_argument(
        "--vcf",
        help="VCF the VEP annotations were called from; when given, each report "
        "line carries the variant's per-allele read depths (e.g. 'T:0 C:562')",
    )
    parser.add_argument(
        "--suppress",
        help="comma-delimited consequence type(s) whose rows are dropped",
    )
    parser.add_argument(
        "--feature_qualifier",
        default="product",
        help="comma-delimited attribute key(s), matched case-insensitively, read "
        "off the CDS descendant(s) to name the product",
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
    suppress = set(split_qualifiers(args.suppress))
    qualifiers = split_qualifiers(args.feature_qualifier)
    depth_index = build_depth_index(args.vcf) if args.vcf else None

    # query labels: explicit --query_genes take precedence over BED names, matching
    # how extract_variants selected the variants being reported on
    if args.query_genes:
        query_list = ordered_query_genes(args.query_genes)
    elif args.bedfile:
        query_list = sorted(extract_queries_from_bed(args.bedfile))
    else:
        query_list = []

    lines = report_variants(
        args.vep_tsv,
        features,
        suppress,
        args.feature_type,
        qualifiers,
        depth_index,
        query_list,
        args.exact_match,
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
