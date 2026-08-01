"""Query-gene identifier resolution shared by the theiagene commands."""

from collections import defaultdict

from theiagene.lib.feature import FeatureCol


def exact_check(query_set: set, id: str) -> bool:
    """Return True or False for an exact match"""
    return id in query_set


def substring_check(query_set: set, id: str) -> bool:
    """Return True or False for a substring match"""
    return any(query in id for query in query_set)


def extract_queries_from_bed(bedfile: str) -> set:
    """Extract query regions from BED"""
    with open(bedfile, "r") as raw:
        return set(x.split()[3] for x in raw)


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
    explicit gene name over the raw feature id."""
    for key in ("gene", "Name", "product"):
        value = feature.attributes.get(key)
        if value:
            return value
    return feature.fid


def gff_query_ranges(
    features: FeatureCol,
    query_list: list,
    group_by: str,
    feature_type: str,
    feature_qualifiers: list,
    exact_match: bool,
) -> dict:
    """Flatten the ``feature_type`` coordinates of the selected query genes to
    ``{<CONTIG>: [(START, END, LABEL), ...]}`` (0-based, half-open).

    The query units are the ``group_by`` features (e.g. each RNA). When
    ``query_list`` is non-empty a unit is kept only if one of its identifiers
    matches a query term, and the matched term becomes the annotation label;
    otherwise every unit is kept and labelled by its own gene name."""
    contig2ranges = defaultdict(list)
    for feature in features[group_by]:
        if query_list:
            label = match_query(
                query_list, feature_identifiers(feature, feature_qualifiers), exact_match
            )
            if label is None:
                continue
        else:
            label = feature_label(feature)
        contig = feature.seqid
        # only the requested-type subfeatures beneath this query unit (e.g. CDS)
        subfeatures = FeatureCol(list(iter_descendants(feature)), group=False)
        for subfeature in subfeatures[feature_type]:
            start, end = subfeature.start, subfeature.end
            if end <= start:
                raise ValueError(
                    f"Invalid region for query '{label}' on contig '{contig}': "
                    f"start ({start}) must be < end ({end})"
                )
            contig2ranges[contig].append((start, end, label))
    return contig2ranges


def bed_query_ranges(bedfile: str, query_set: set) -> dict:
    """Flatten a BED file to ``{<CONTIG>: [(START, END, NAME), ...]}`` (0-based,
    half-open), keeping only rows whose name is in ``query_set`` (an empty set
    keeps every row).

    A BED region is used directly as a query coordinate; the name column (col 4)
    is both the filter key and the annotation label."""
    contig2ranges = defaultdict(list)
    with open(bedfile) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            data = line.split()
            if not data:
                continue
            contig, start, end, name = data[0], int(data[1]), int(data[2]), data[3]
            if query_set and name not in query_set:
                continue
            if end <= start:
                raise ValueError(
                    f"Invalid region for query '{name}' on contig '{contig}': "
                    f"start ({start}) must be < end ({end})"
                )
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
