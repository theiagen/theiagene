"""Render SnpEff variant annotations into a gene-labelled report TSV.

Given a SnpEff-annotated VCF and the reference GFF used to build its database,
this command keeps the ``ANN`` entries whose consequence is not suppressed,
resolves each entry's variant back to a CDS product name through the GFF feature
hierarchy, and leads the HGVS changes with a single gene label, e.g.::

    ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)

The label is the ``--query_genes`` term that matches the entry's feature -- the
name the caller asked about, not the full product it resolved to -- and the
product follows in quotes so a name carrying commas survives being joined into
a comma-delimited report. Without ``--query_genes`` (or when none match) the
label falls back to the normalized product name. The trailing per-allele read
depths are read off the annotated record's own sample (``AD``, else
``RO``/``AO``) and are omitted when it carries neither.

An entry is tied back to the annotation by the variant's coordinates; its
``Feature_ID`` is consulted only to break a tie between several annotation
units overlapping one variant.

Each kept entry becomes one record of a TSV (stdout, or ``--output``) whose
``#``-prefixed header is :data:`COLUMNS`: the gene label, the HGVS.c and HGVS.p
strings prefixed with their transcript and protein identifiers, their
abbreviated nucleotide (``428A>G``) and amino acid (``K143R``) changes, and the
report line above. A field the entry has no HGVS string for is ``NA``, e.g.::

    #GENE  HGVSc               HGVSp                 NT      AA     REPORT
    ERG11  rna-x:c.428A>G      prot-x:p.Lys143Arg    428A>G  K143R  ERG11: "..." (...)

HGVS strings and consequence terms are SnpEff's own, so they inherit its
notation (see the README's deviations from HGVS)."""

import re
import sys
import logging
import argparse

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

# the '|'-delimited subfields of one SnpEff ANN entry, in specification order
_ANN_FIELDS = (
    "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
    "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c",
    "HGVS.p", "cDNA.pos", "CDS.pos", "AA.pos", "Distance", "Messages",
)

# SnpEff flags a transcript whose CDS lacks a start codon or a whole number of
# codons; its translation is unreliable, so such a transcript is treated as
# noncoding and its protein change dropped
_PARTIAL_TRANSCRIPT = {"WARNING_TRANSCRIPT_INCOMPLETE", "WARNING_TRANSCRIPT_NO_START_CODON"}

# the attribute keys a transcript can be identified by; tried in turn when
# several units overlap one variant
_UNIT_ID_KEYS = ("ID", "Name", "transcript_id")


def iter_annotations(vcf: str):
    """Yield ``(record, entry)`` for every SnpEff ``ANN`` entry of a VCF, where
    ``entry`` maps :data:`_ANN_FIELDS` to that entry's subfields.

    Records carrying no ``ANN`` (e.g. intergenic variants under
    ``-no-intergenic``) yield nothing."""
    handle = import_vcf(vcf)
    try:
        if "ANN" not in handle.header.info:
            raise ValueError(f"{vcf} carries no SnpEff ANN annotations")
        for record in handle:
            for raw in record.info.get("ANN") or ():
                yield record, dict(zip(_ANN_FIELDS, raw.split("|")))
    finally:
        handle.close()


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


def _depth_suffix(record, allele: str):
    """Render ``"{ref}:{ref_depth} {alt}:{alt_depth}"`` for ``allele`` of a VCF
    record, or None when the record carries no per-allele depths or ``allele``
    is not one of its alts."""
    sample = next(iter(record.samples.values()), None)
    depths = _allele_depths(sample) if sample is not None else None
    alts = record.alts or ()
    if depths is None or allele not in alts:
        return None
    ref_depth, alt_depths = depths
    return f"{record.ref}:{ref_depth} {allele}:{alt_depths[alts.index(allele)]}"


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


# HGVS three-letter amino acid codes -> their one-letter abbreviations
_AA_CODES = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Xaa": "X", "Ter": "*",
}
_THREE_LETTER = re.compile("|".join(_AA_CODES))


def _abbreviate(hgvs: str) -> str:
    """Collapse an HGVS change into its abbreviated, prefix-free form.

    Nucleotide changes keep their notation (``c.428A>G`` -> ``428A>G``,
    ``c.383delA`` -> ``383delA``), whereas amino acid changes have their
    three-letter residues replaced by one-letter codes (``p.Lys143Arg`` ->
    ``K143R``, ``p.Arg458*`` -> ``R458*``, ``p.Lys128fs`` -> ``K128fs``)."""
    body = hgvs.split(".", 1)[-1]
    if hgvs.startswith("p."):
        return _THREE_LETTER.sub(lambda match: _AA_CODES[match.group()], body)
    return body


def report_line(label: str, product: str, consequence: str, hgvsc: str = None,
                hgvsp: str = None, depths: str = None) -> str:
    """Build the report line for a single kept annotation.

    The line leads with ``label`` -- the query-gene term the variant was
    extracted under (e.g. ``FKS1``) -- followed by the quoted CDS product; the
    consequence and whichever HGVS changes are present follow in parentheses.
    Quoting the product keeps a name carrying commas ('1,3-beta-glucan synthase
    ...') intact once the caller joins the lines with ','. ``depths`` (see
    :func:`_depth_suffix`) is appended after a ``;`` when given (e.g. ``...
    p.Lys143Arg; T:0 C:562``)."""
    body = " ".join(piece for piece in (consequence, hgvsc, hgvsp) if piece)
    if depths:
        body += f"; {depths}"
    return f'{label}: "{product}" ({body})'


def _select_unit(units: list, feature_id: str):
    """Pick the single annotation unit an annotation is about, or None.

    An annotation's HGVS strings are computed against exactly one transcript,
    and SnpEff has already emitted a separate entry per transcript, so an entry
    must resolve to one unit rather than fan out over every unit the variant
    overlaps. A sole overlapping unit is taken as-is -- the coordinate has
    already done the work and there is nothing to disambiguate. Only where
    several overlap does the entry's ``Feature_ID`` break the tie, compared
    against the unit's own identifiers (never its parent's or its CDS product's,
    which sibling transcripts of one gene share and which would therefore match
    both).

    None is returned when nothing overlapped, when an ambiguous overlap carries
    no ``Feature_ID`` to arbitrate it, and when none of the units answer to the
    one it carries -- all of which leave the caller to skip the entry rather
    than attribute it to whichever unit happened to sort first."""
    if len(units) < 2:
        return units[0] if units else None
    if not feature_id:
        return None
    for unit in units:
        # `fid` alongside the attributes: group_features may have renamed a
        # colliding ID, leaving the one SnpEff saw only in `attributes`
        identifiers = [unit.fid] + [unit.attributes.get(key) for key in _UNIT_ID_KEYS]
        if any(identifier == feature_id for identifier in identifiers):
            return unit
    return None


# report TSV columns: the gene label, the identifier-prefixed HGVS strings,
# their abbreviations and the formatted report line
COLUMNS = ("GENE", "HGVSc", "HGVSp", "NT", "AA", "REPORT")


def report_variants(
    vcf: str,
    features: FeatureCol,
    suppress: set,
    feature_type: str,
    qualifiers: list,
    query_list: list = None,
    exact_match: bool = False,
) -> list:
    """Turn a SnpEff-annotated VCF into query-labelled report records, one per
    kept ``ANN`` entry, whose fields follow :data:`COLUMNS`.

    An entry is dropped when any of its consequence terms is suppressed, when it
    is not annotated against a transcript (e.g. ``intergenic_region``), when it
    carries neither an HGVS.c nor an HGVS.p string, or when its variant cannot
    be resolved to a CDS product in ``features`` -- either because no annotation
    unit overlaps it, because several do and none answers to the entry's
    ``Feature_ID`` (see :func:`_select_unit`), or because the resolved unit
    carries no qualifier. Each kept record is labelled with the ``query_list``
    term that matched the entry's unit, falling back to the product-derived
    label when no query matched.

    HGVS.c is prefixed with the transcript's ``Feature_ID`` and HGVS.p with its
    CDS ``protein_id`` (else the ``Feature_ID``); NT/AA are their abbreviations
    (see :func:`_abbreviate`), and a field the entry has no HGVS string for is
    ``NA``. A transcript SnpEff flags as partial (:data:`_PARTIAL_TRANSCRIPT`) is
    treated as noncoding, so its HGVS.p is dropped."""
    records = []
    for record, entry in iter_annotations(vcf):
        if any(term in suppress for term in entry["Annotation"].split("&")):
            continue
        if entry["Feature_Type"] != "transcript":
            continue
        hgvsc = entry["HGVS.c"]
        hgvsp = entry["HGVS.p"]
        if _PARTIAL_TRANSCRIPT & set(entry["Messages"].split("&")):
            hgvsp = ""
        if not hgvsc and not hgvsp:
            continue
        variant = f"{record.contig}:{record.pos} {record.ref}>{entry['Allele']}"
        # units are taken from one level: the transcripts SnpEff annotates
        # against, else the genes of an annotation carrying no transcript level
        # (gene -> CDS), else bare feature_type records owning no gene at all.
        # Taking one level keeps a locus from counting once per tier of its own
        # hierarchy
        hits = features.index(record.contig, record.start, record.stop)
        units = hits.rnas or hits.genes or hits[feature_type]
        feature = _select_unit(units, entry["Feature_ID"])
        qualifier_hit = _type_qualifier(feature, feature_type, qualifiers)
        if qualifier_hit is None:
            logger.warning(
                f"no {feature_type} {qualifiers} qualifier resolved for "
                f"{entry['Feature_ID']}; skipping variant {variant}"
            )
            continue
        label = _query_label(feature, query_list, qualifiers, exact_match)
        if query_list and label is None:
            logger.warning(
                f"no query gene matched {entry['Feature_ID']}; labelling variant "
                f"{variant} by its product instead"
            )
        label = label or normalize_name(qualifier_hit)
        protein_id = _type_qualifier(feature, feature_type, ["protein_id"]) or entry["Feature_ID"]
        records.append([
            label,
            f"{entry['Feature_ID']}:{hgvsc}" if hgvsc else "NA",
            f"{protein_id}:{hgvsp}" if hgvsp else "NA",
            _abbreviate(hgvsc) if hgvsc else "NA",
            _abbreviate(hgvsp) if hgvsp else "NA",
            report_line(label, qualifier_hit, entry["Annotation"], hgvsc, hgvsp,
                        _depth_suffix(record, entry["Allele"])),
        ])
    return records


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register the report_variants arguments on ``parser``"""
    parser.add_argument(
        "--vcf",
        required=True,
        help="SnpEff-annotated VCF; per-allele read depths (e.g. 'T:0 C:562') "
        "are read off its sample",
    )
    parser.add_argument("--reference_gff", required=True)
    parser.add_argument(
        "--query_genes",
        help="comma-delimited query gene name(s)"
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
        "--suppress",
        help="comma-delimited SnpEff consequence term(s) whose annotations are dropped",
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
        help="write the report TSV here instead of stdout",
    )
    return parser


def run_cli(args: argparse.Namespace) -> int:
    """Render the SnpEff annotations into a product-named report TSV"""
    features = assimilate_gff(args.reference_gff)
    suppress = set(split_qualifiers(args.suppress))
    qualifiers = split_qualifiers(args.feature_qualifier)

    # query labels: explicit --query_genes take precedence over BED names, matching
    # how extract_variants selected the variants being reported on
    if args.query_genes:
        query_list = ordered_query_genes(args.query_genes)
    elif args.bedfile:
        query_list = sorted(extract_queries_from_bed(args.bedfile))
    else:
        query_list = []

    records = report_variants(
        args.vcf,
        features,
        suppress,
        args.feature_type,
        qualifiers,
        query_list,
        args.exact_match,
    )

    tsv = "\n".join(["#" + "\t".join(COLUMNS)] + ["\t".join(record) for record in records])
    if args.output:
        with open(args.output, "w") as handle:
            handle.write(tsv + "\n")
        logger.debug(f"Wrote {len(records)} variant record(s) to {args.output}")
    else:
        print(tsv)
        logger.debug(f"Reported {len(records)} variant(s)")

    return 0


def main(argv=None) -> int:
    """Standalone entrypoint (``python -m theiagene.report_variants``)"""
    parser = argparse.ArgumentParser(
        description="Extracts variant annotations from a SnpEff-annotated VCF"
    )
    add_arguments(parser)
    args = parser.parse_args(argv)
    configure_logging()
    return run_cli(args)


if __name__ == "__main__":
    sys.exit(main())
