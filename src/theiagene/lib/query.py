"""Query-gene identifier resolution shared by the theiagene commands.

Both ``gene_coverage`` and ``variant_annotation`` accept a set of query genes
and must decide whether a feature identifier (product name, locus tag, ...)
matches one of them.  This module collects the identifier-matching and
normalization helpers used by both."""


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
    'lanosterol.14-alpha.demethylase') match spaced GBFF products
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
