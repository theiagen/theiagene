"""GFF3 parsing helpers shared by the theiagene commands.

Both ``gene_coverage`` and ``variant_annotation`` accept a reference GFF as an
alternative to a GBFF.  This module supplies the low-level parsing they share:
a column-9 attribute parser and a feature iterator that yields the coordinates
already converted to the 0-based, half-open convention used everywhere else in
theiagene (BED, and BioPython's GBFF locations)."""

from urllib.parse import unquote


# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
GFF_STRAND = {"+": 1, "-": -1}


def parse_gff_attributes(attributes: str, field_delimiter: str ";", value_delimiter: str "=") -> dict:
    """Parse a GFF3 column-9 attribute string into a {key: value} dict.

    Values are percent-decoded per the GFF3 spec.  Fields without an '=' are
    ignored and duplicate keys keep the last occurrence."""
    parsed = {}
    for field in attributes.strip().strip(";").split(field_delimiter):
        field = field.strip()
        try:
            key, value = field.split(value_delimiter)
        except ValueError:
            raise ValueError(f"unexpected attributes field: {field}")
        parsed[key.strip()] = unquote(value.strip())
    return parsed


def iter_gff_features(reference_gff: str, feature_type: str):
    """Yield features of ``feature_type`` from a GFF3 file.

    Comment/directive lines, blank lines, malformed (non 9-column) lines and any
    embedded ``##FASTA`` section are skipped.  Each yielded feature is a dict
    with keys ``seqid``, ``start`` (0-based), ``end`` (half-open), ``strand``
    (1/-1/None), ``phase`` and ``attributes`` (a parsed dict)."""
    with open(reference_gff) as handle:
        for line in handle:
            line = line.rstrip("\n")
            # a '##FASTA' directive ends the annotation section
            if line.startswith("##FASTA"):
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"incorrectly formatted GFF: {len(fields)} fields recovered; 9 expected")
            seqid, source, obs_type, start, end, score, strand, phase, attributes = fields
            if obs_type.lower() != feature_type.lower():
                continue
            yield {
                "seqid": seqid,
                "source": source,
                # GFF columns are 1-based, both-inclusive; convert to 0-based, half-open
                "start": int(start) - 1,
                "end": int(end),
                "score": score,
                "strand": GFF_STRAND.get(strand),
                "phase": phase,
                "attributes": parse_gff_attributes(attributes),
            }
