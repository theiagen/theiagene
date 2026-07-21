"""Shared reference-parsing layer for the theiagene commands.

Both ``gene_coverage`` and ``variant_annotation`` read query-gene coordinates
from the same reference formats (GenBank/GFF3, plus BED for coverage).  This
module holds the reading, multi-exon grouping and identifier-matching they share,
so each command keeps only its own object-construction step:

- :func:`iter_gbff_raw` / :func:`iter_gff_raw` yield a neutral :class:`RawGene`
  per candidate feature (pre-match, no naming) -- ``gene_coverage`` turns these
  into :class:`~theiagene.lib.gene_model.Gene` objects, ``variant_annotation``
  into :class:`~theiagene.lib.gene_model.GeneModel` objects.
- :func:`match_identifiers` reconciles the two commands' divergent matching
  behind ``normalize`` / ``exact_match`` / ``id_qualifiers`` options.
- :func:`resolve_contig` applies the contig-membership check (after a gene has
  matched, so unmatched genes on foreign contigs never raise).
- :func:`parse_bed_genes` handles the coverage-only BED format directly.
"""

from collections import namedtuple

from Bio import SeqIO

from theiagene.lib.gene_model import Gene
from theiagene.lib.gff import iter_gff_features
from theiagene.lib.query import exact_check, substring_check, match_query


# a neutral, pre-match feature: coordinates + metadata, no query naming applied.
# ``qualifiers`` is a {name: [values]} dict (uniform across GBFF and GFF);
# ``contig_seq`` is the reference sequence when available, else None.
RawGene = namedtuple(
    "RawGene", "contig_candidates strand qualifiers parts contig_seq"
)


# GFF3 attribute keys used to coalesce multi-segment (multi-exon) CDS lines into
# one gene, in preference order; all segments of a CDS share these (Parent/ID
# especially)
_GFF_GROUP_KEYS = ("Parent", "ID", "locus_tag", "gene", "protein_id")


def iter_gbff_raw(reference_gbff: str, feature_type: str):
    """Yield a :class:`RawGene` for every ``feature_type`` feature in a GBFF"""
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_seq = str(record.seq)
            for feature in record.features:
                if feature.type.lower() != feature_type.lower():
                    continue
                # GenBanks are 1-based; BioPython adjusts to 0-based half-open
                parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
                yield RawGene(
                    contig_candidates=[record.id, record.name],
                    strand=feature.location.strand,
                    qualifiers=dict(feature.qualifiers),
                    parts=parts,
                    contig_seq=contig_seq,
                )


def iter_gff_raw(reference_gff: str, feature_type: str, fa_dict: dict = None):
    """Yield a :class:`RawGene` per gene in a GFF3, coalescing multi-exon CDS.

    A multi-exon CDS is spread over several GFF lines sharing an ``ID``/``Parent``;
    those are collapsed into one :class:`RawGene` carrying every segment, in
    first-seen order.  ``qualifiers`` accumulates the attribute values across all
    of a gene's lines.  When ``fa_dict`` (a ``{seqid: SeqRecord}`` mapping) is
    given, ``contig_seq`` is filled from it (``None`` when the seqid is absent, so
    the caller can raise with a format-appropriate message); otherwise it is
    ``None`` (coordinate-only use)."""
    groups = {}
    order = []
    for feature in iter_gff_features(reference_gff, feature_type):
        attrs = feature["attributes"]
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
                "qualifiers": {},
                "parts": [],
            }
            groups[group_id] = group
            order.append(group_id)
        group["parts"].append((feature["start"], feature["end"]))
        # accumulate attribute values across the gene's lines (uniform list form)
        for name, value in attrs.items():
            values = group["qualifiers"].setdefault(name, [])
            if value not in values:
                values.append(value)

    for group_id in order:
        group = groups[group_id]
        seqid = group["seqid"]
        yield RawGene(
            contig_candidates=[seqid],
            strand=group["strand"],
            qualifiers=group["qualifiers"],
            parts=group["parts"],
            contig_seq=(
                str(fa_dict[seqid].seq)
                if fa_dict is not None and seqid in fa_dict
                else None
            ),
        )


def match_identifiers(
    qualifiers: dict,
    query_list,
    id_qualifiers,
    exact_match: bool = False,
    normalize: bool = True,
):
    """Match a feature's qualifiers against the query list.

    Returns ``(matched, identifiers)`` where ``identifiers`` is every value
    collected across ``id_qualifiers`` (in order).  When ``normalize`` is True
    (variant_annotation) matching is normalization-aware via
    :func:`theiagene.lib.query.match_query` and ``matched`` is the query term that
    matched; when False (gene_coverage) a raw exact/substring test is applied per
    identifier and ``matched`` is the matching identifier itself.  ``matched`` is
    ``None`` when nothing matches."""
    identifiers = []
    for qualifier in id_qualifiers:
        identifiers.extend(qualifiers.get(qualifier, []))
    if not identifiers:
        return None, identifiers
    if normalize:
        return match_query(query_list, identifiers, exact_match), identifiers
    query_set = set(query_list)
    check = exact_check if exact_match else substring_check
    matched = next((ident for ident in identifiers if check(query_set, ident)), None)
    return matched, identifiers


def resolve_contig(candidates, contig_names, require: bool, source_label: str) -> str:
    """Return the first candidate present in ``contig_names``.

    Falls back to the first candidate when ``require`` is False; raises
    ``KeyError`` when required and none of the candidates are present."""
    for candidate in candidates:
        if candidate in contig_names:
            return candidate
    if require:
        raise KeyError(f"{' and '.join(candidates)} not in {source_label}")
    return candidates[0]


def parse_bed_genes(
    bedfile: str,
    query_list,
    exact_match: bool,
    contig_names,
    require: bool = True,
    source_label: str = "BAM",
) -> list:
    """Parse a BED file into :class:`Gene` objects (coverage-only format).

    Rows sharing a name on the same contig accumulate as multiple parts, matching
    the multi-segment handling of the GBFF/GFF parsers."""
    query_set = set(query_list)
    check = exact_check if exact_match else substring_check
    genes = {}
    order = []
    with open(bedfile) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            data = line.split()
            if not data:
                continue
            name = data[3]
            if not check(query_set, name):
                continue
            # BED files are already 0-based, half-open
            contig = resolve_contig([data[0]], contig_names, require, source_label)
            key = (contig, name)
            gene = genes.get(key)
            if gene is None:
                gene = Gene(gene_id=name, contig=contig)
                genes[key] = gene
                order.append(key)
            gene.add_part(int(data[1]), int(data[2]))
    return [genes[key] for key in order]
