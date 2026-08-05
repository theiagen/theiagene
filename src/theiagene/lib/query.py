"""Query-gene identifier resolution shared by the theiagene commands."""

from collections import defaultdict

from theiagene.lib.feature import FeatureCol
from theiagene.lib.parsers import iter_bed_rows


def exact_check(query_set: set, id: str) -> bool:
    """Return True or False for an exact match"""
    return id in query_set


def substring_check(query_set: set, id: str) -> bool:
    """Return True or False for a substring match"""
    return any(query in id for query in query_set)


def extract_queries_from_bed(bedfile: str) -> set:
    """Extract query regions from BED (the name column of each data row)"""
    return {data[3] for data in iter_bed_rows(bedfile)}


def normalize_name(name: str) -> str:
    """Collapse whitespace and VCF-reserved characters to '.' for a stable id.

    e.g. 'lanosterol 14-alpha demethylase' -> 'lanosterol.14-alpha.demethylase'"""
    out = name
    for char in (" ", "\t", ";", "=", ","):
        out = out.replace(char, ".")
    # collapse runs of dots that can arise from adjacent replaced characters
    while ".." in out:
        out = out.replace("..", ".")
    return out.strip(".")


def sanitize_info_value(value: str) -> str:
    """Sanitize a string for use in a VCF INFO field (no whitespace or reserved characters)"""
    for char in (" ", "\t", ";", "=", ","):
        value = value.replace(char, "_")
    return value


def match_query(query_list, identifiers, exact_match: bool):
    """Return the first query term (in order) matching any candidate identifier.

    Matching is normalization-aware, so dotted query names (e.g.
    'lanosterol.14-alpha.demethylase') match spaced products
    ('lanosterol 14-alpha demethylase') and vice versa."""
    for query in query_list:
        nq = normalize_name(query)
        for ident in identifiers:
            ni = normalize_name(ident)
            if exact_match:
                if query == ident or nq == ni:
                    return query
            elif query in ident or nq in ni:
                return query
    return None


def ordered_query_genes(query_genes_arg) -> list:
    """Flatten the --query_genes argument into an ordered, de-duplicated list"""
    ordered = []
    seen = set()
    for chunk in query_genes_arg or []:
        for gene in chunk.split(","):
            gene = gene.strip()
            if gene and gene not in seen:
                seen.add(gene)
                ordered.append(gene)
    return ordered


def split_qualifiers(raw) -> list:
    """Split a comma-/space-delimited qualifier string into individual keys,
    dropping empty tokens. A None/empty input yields an empty list."""
    if not raw:
        return []
    return raw.replace(",", " ").split()


def iter_descendants(feature):
    """Yield every descendant of ``feature``, depth-first (children of children
    included), so a query unit's CDS are reachable no matter how deep the
    gene -> RNA -> CDS hierarchy runs."""
    for descendant in feature.descendants:
        yield descendant
        yield from iter_descendants(descendant)


def feature_identifiers(feature, qualifiers) -> list:
    """Collect the candidate name strings a query term may match against.

    A gene can be named on the query feature itself, on its parent (the gene
    record usually carries ``gene``/``Name``) or on its CDS descendants (which
    usually carry ``product``); every such identifier is gathered so a query
    matches regardless of which record holds the name. The attribute keys to
    read are supplied by ``qualifiers`` and matched case-insensitively."""
    identifiers = []
    related = [feature]
    if feature.parent is not None:
        related.append(feature.parent)
    related.extend(iter_descendants(feature))
    wanted = {qualifier.lower() for qualifier in qualifiers}
    for related_feature in related:
        if related_feature.fid:
            identifiers.append(related_feature.fid)
        for key, value in related_feature.attributes.items():
            if value and key.lower() in wanted:
                identifiers.append(value)
    return identifiers


def feature_label(feature) -> str:
    """Return a stable human-readable name for a query feature, preferring an
    explicit gene name over the raw feature id.

    The name may live on the query feature itself, on its parent (the gene
    record usually carries ``gene``/``Name``) or on a CDS descendant (which
    usually carries ``product``); the same related features consulted by
    :func:`feature_identifiers` are searched so distinct genes resolve to
    distinct labels."""
    related = [feature]
    if feature.parent is not None:
        related.append(feature.parent)
    related.extend(iter_descendants(feature))
    for key in ("gene", "Name", "product"):
        for related_feature in related:
            value = related_feature.attributes.get(key)
            if value:
                return value
    return feature.fid


def validate_region(start: int, end: int, label: str, contig: str, contig_len: int = None) -> None:
    """Raise ``ValueError`` when a query region's coordinates are invalid.

    Query ranges are 0-based, half-open, so ``start`` must be < ``end``; this is
    checked wherever a range is built (GFF and BED) and again at quantification.
    When a ``contig_len`` is supplied (the BAM stage, which knows the reference
    length) the region must also fall within ``[0, contig_len]``. The message
    names the query label and contig so a bad region is traceable back to the
    GFF record or BED row it came from."""
    prefix = f"Invalid region for query '{label}' on contig '{contig}': "
    if end <= start:
        raise ValueError(prefix + f"start ({start}) must be < end ({end})")
    if contig_len is not None:
        if start < 0:
            raise ValueError(prefix + f"start ({start}) must be >= 0")
        if end > contig_len:
            raise ValueError(prefix + f"end ({end}) exceeds contig length ({contig_len})")


def _grouped_query_ranges(
    features: FeatureCol,
    query_list: list,
    group_by: str,
    feature_type: str,
    feature_qualifiers: list,
    exact_match: bool,
) -> dict:
    """Flatten the ``feature_type`` coordinates of the query units selected at
    the ``group_by`` level to ``{<CONTIG>: [(START, END, LABEL), ...]}``
    (0-based, half-open); see :func:`gff_query_ranges`."""
    contig2ranges = defaultdict(list)
    for feature in features[group_by]:
        if query_list:
            matched = match_query(
                query_list, feature_identifiers(feature, feature_qualifiers), exact_match
            )
            if matched is None:
                continue
        label = feature_label(feature)
        contig = feature.seqid
        # the requested-type features at or beneath this query unit (e.g. CDS);
        # the unit itself is included so that grouping directly on the
        # feature_type (the fallback level) still resolves its own coordinates
        subfeatures = FeatureCol([feature] + list(iter_descendants(feature)), group=False)
        for subfeature in subfeatures[feature_type]:
            start, end = subfeature.start, subfeature.end
            validate_region(start, end, label, contig)
            contig2ranges[contig].append((start, end, label))
    return contig2ranges


def gff_query_ranges(
    features: FeatureCol,
    query_list: list,
    feature_type: str,
    feature_qualifiers: list,
    exact_match: bool,
) -> dict:
    """Flatten the ``feature_type`` coordinates of the selected query genes to
    ``{<CONTIG>: [(START, END, LABEL), ...]}`` (0-based, half-open).

    The query units are resolved by trying the ``gene`` level first and falling
    back to the ``feature_type`` level (e.g. ``CDS``) when grouping by gene
    resolves nothing -- so a bacterial annotation whose ``CDS`` is a direct child
    of ``gene`` (no intervening ``mRNA``) still yields coordinates, as does an
    annotation carrying only bare ``feature_type`` records with no gene parent.

    When ``query_list`` is non-empty a unit is kept only if one of its
    identifiers matches a query term; the query term is used only to filter, and
    every kept unit is labelled by its own resolved gene name so distinct genes
    remain distinct rows (a query matching two paralogs yields two labels)."""
    # fallback order: the gene that owns the feature_type, then the raw
    # feature_type itself for annotations that carry no gene record
    for group_by in ("gene", feature_type):
        contig2ranges = _grouped_query_ranges(
            features, query_list, group_by, feature_type, feature_qualifiers, exact_match
        )
        if contig2ranges:
            return contig2ranges
    return contig2ranges


def bed_query_ranges(bedfile: str, query_set: set) -> dict:
    """Flatten a BED file to ``{<CONTIG>: [(START, END, NAME), ...]}`` (0-based,
    half-open), keeping only rows whose name is in ``query_set`` (an empty set
    keeps every row).

    A BED region is used directly as a query coordinate; the name column (col 4)
    is both the filter key and the annotation label."""
    contig2ranges = defaultdict(list)
    for data in iter_bed_rows(bedfile):
        contig, start, end, name = data[0], int(data[1]), int(data[2]), data[3]
        if query_set and name not in query_set:
            continue
        validate_region(start, end, name, contig)
        contig2ranges[contig].append((start, end, name))
    return contig2ranges


def collapse_to_single_contig(contig2ranges: dict, contig: str) -> dict:
    """Re-file every range under a single ``contig`` (the ambiguous-contig case,
    where the sample was mapped to a reference whose contig name differs from
    the coordinate source)."""
    collapsed = defaultdict(list)
    for ranges in contig2ranges.values():
        collapsed[contig].extend(ranges)
    return collapsed
