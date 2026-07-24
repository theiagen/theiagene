"""Feature/gene data model shared by the theiagene commands."""

from io import StringIO
from collection import defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
_STRAND = {"+": 1, "-": -1}
# strand integer -> GFF strand column (anything else, e.g. None, serializes to '.')
_STRAND_REVERSE = {1: "+", -1: "-"}
# the placeholder characters GFF3 uses for an undefined numeric/strand column
_UNDEFINED = {".", "?", "", None}

# GFF3 column-9 reserved characters and their percent-encodings; '%' is listed
# first so the escapes introduced below are not themselves re-encoded
_GFF_ATTR_ESCAPES = (
    ("%", "%25"), (";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C"),
    ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D"),
)


def _gff_escape(value: str) -> str:
    """Percent-encode the GFF3 attribute-reserved characters in ``value``"""
    for char, code in _GFF_ATTR_ESCAPES:
        value = value.replace(char, code)
    return value


def _format_gff_attributes(attributes: dict, field_delimiter: str = ";", value_delimiter: str = "=") -> str:
    """Serialize a ``{key: value}`` attribute dict to a GFF3 column-9 string.

    Reserved characters in keys and values are percent-encoded; an empty dict
    yields ``.`` (the GFF3 "no attributes" placeholder)."""
    if not attributes:
        return "."
    return field_delimiter.join(
        f"{_gff_escape(str(key))}{value_delimiter}{_gff_escape(str(value))}"
        for key, value in attributes.items()
    )


def _as_int(value, name: str):
    """Coerce a GFF numeric column to ``int``, mapping the undefined placeholder
    ('.') to ``None``"""
    if value in _UNDEFINED:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be an integer, got {value!r}")


def _as_strand(value):
    """Coerce a GFF strand column ('+'/'-') or a signed integer to 1/-1, with the
    undefined placeholder ('.'/'?') mapping to ``None``"""
    if value in _UNDEFINED:
        return None
    if value in _STRAND:
        return _STRAND[value]
    if value in (1, -1):
        return value
    raise ValueError(f"strand must be one of '+'/'-'/1/-1, got {value!r}")


class Feature():
    """A single cross-annotation feature
    `start` and `end` are sorted by smallest to largest and are expected to be 0-based half-open upon input"""

    def __init__(self, fid=None, pid=None, seqid=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=None,
                 attributes: dict = None, descendants: list = None, parent: Feature = None, sequence: str = "", origin_sequence: str = ""
                 ):
        self.id = fid
        # parent ID
        self.pid = pid
        self.seqid = seqid
        self.source = source
        self.type = type
        put_start = _as_int(start, "start")
        put_end = _as_int(end, "end")
        # enforce start < end
        self.start, self.end = sorted((put_start, put_end))
        self.score = _as_int(score, "score")
        self.strand = _as_strand(strand)
        self.phase = _as_int(phase, "phase")
        self.parent = parent

        # sentinel defaults: a shared mutable default would let every Feature
        # alias one dict/list, collapsing the hierarchy built by group_features
        self.attributes = {} if attributes is None else attributes
        self.descendants = [] if descendants is None else descendants

        # derive sequence from the full contiguous sequence provided
        if origin_sequence:
            # 0-BASED, HALF-OPEN
            sequence = origin_sequence[self.start:self.end]

        # demand that the provided sequence abide by the coordinates
        if sequence:
            if len(sequence) < self.end - self.start:
                raise ValueError(f"sequence length deviates from coordinate length")
            else:
                self.sequence = sequence
        else:
            self.sequence = ""


    def _to_gff_line(self) -> str:
        """Serialize this feature (without its descendants) to one GFF3 line.

        Coordinates are converted from the internal 0-based, half-open
        representation back to GFF3's 1-based, both-inclusive columns; an
        undefined ``start``/``end``/``score``/``strand``/``phase`` renders as the
        '.' placeholder."""
        start = "." if self.start is None else str(self.start + 1)
        end = "." if self.end is None else str(self.end)
        score = "." if self.score is None else str(self.score)
        strand = _STRAND_REVERSE.get(self.strand, ".")
        phase = "." if self.phase is None else str(self.phase)
        id_dict = {"ID": self.fid}
        if self.pid:
            id_dict["Parent"] = self.pid
        return "\t".join((
            self.seqid, self.source, self.type, start, end,
            score, strand, phase, _format_gff_attributes({**id_dict, **self.attributes)),
        ))
    
    def to_gff(self) -> str:
        """Serialize this feature and all of its descendants to a GFF3 string.

        Lines are emitted depth-first, each parent before its children, joined
        by newlines (no trailing newline)."""
        lines = [self._to_gff_line()]
        for descendant in self.descendants:
            lines.append(descendant.to_gff())
        return "\n".join(lines)


def _attr_get(attributes: dict, keys):
    """Return the value of the first present attribute key (accepts case
    variants, e.g. ``ID``/``id`` or ``Parent``/``parent``), or ``None``"""
    for key in keys:
        if key in attributes:
            return attributes[key]
    return None


def group_features(features: list, parent_ids: list = ["Parent", "parent"], ids: list = ["ID", "Id", "id"]):
    """Hierarchically group features based on Parent <-> ID relationships"""
    id2feature = {}
    for feature in features:
        fid = _attr_get(feature.attributes, ids)
        if fid in id2feature:
            raise KeyError(f"{fid} is depicted in multiple Features")
        id2feature[fid] = feature

    feature_dict = defaultdict(list)
    for fid, feature in id2feature.items():
        par_id = _attr_get(feature.attributes, parent_ids)
        if par_id:
            # link features together
            id2feature[par_id].descendants.append(feature)
            feature.parent = id2feature[par_id]
        feature_dict[seqid].append(feature)

    return feature_dict