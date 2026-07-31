"""Feature/gene data model shared by the theiagene commands."""

from collections import defaultdict, Counter

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


def _as_int(value, name="value"):
    """Coerce a GFF3 numeric column to an int, or None for an undefined
    ('.'/'?'/''/None) column.

    A value already of int type is returned unchanged. `name` labels the column
    in the error raised when a defined value cannot be parsed as an int."""
    if isinstance(value, int):
        return value
    if value in _UNDEFINED:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} is not an integer: {value!r}")


def _as_strand(value):
    """Coerce a GFF3 strand column ('+'/'-') to a BioPython-style strand integer
    (1/-1), or None for an undefined ('.'/'?'/''/None) strand.

    A value already given as a strand integer (1/-1) is returned unchanged."""
    if value in _STRAND:
        return _STRAND[value]
    if value in (1, -1):
        return value
    if value in _UNDEFINED:
        return None
    raise ValueError(f"strand is not a valid GFF3 strand: {value!r}")


def _unescape_gff_attr(text: str) -> str:
    """Reverse the percent-encoding applied by ``_format_gff_attributes``.

    Escapes are undone in reverse of the encoding order so the '%25' -> '%' step
    runs last and cannot corrupt an escape sequence produced by an earlier
    replacement."""
    for char, escape in reversed(_GFF_ATTR_ESCAPES):
        text = text.replace(escape, char)
    return text


def _parse_gff_attributes(field: str) -> dict:
    """Parse a raw GFF3 column-9 field into an ``{attribute: value}`` dict.

    The field is split on ';' into ``key=value`` pairs; empty segments (from a
    trailing or doubled ';') and bare tokens lacking '=' are skipped. Reserved
    characters percent-encoded by ``_format_gff_attributes`` are decoded in both
    key and value, and a later duplicate key overrides an earlier one. The
    undefined placeholder ('.'/''/None) yields an empty dict."""
    attributes = {}
    if field in _UNDEFINED:
        return attributes
    for pair in field.split(";"):
        if not pair:
            continue
        key, sep, value = pair.partition("=")
        if not sep:
            # a bare token with no '=' carries no value to record
            continue
        attributes[_unescape_gff_attr(key)] = _unescape_gff_attr(value)
    return attributes


def _escape_gff_attr(text: str) -> str:
    """Percent-encode the GFF3 column-9 reserved characters in `text`.

    Replacements run in the order of ``_GFF_ATTR_ESCAPES`` (with '%' first) so a
    '%' introduced by an escape is not itself re-encoded; this is the exact
    inverse of ``_unescape_gff_attr``."""
    for char, escape in _GFF_ATTR_ESCAPES:
        text = text.replace(char, escape)
    return text


def _format_gff_attributes(attributes: dict) -> str:
    """Serialize an ``{attribute: value}`` dict into a GFF3 column-9 field.

    Keys and values have their reserved characters percent-encoded and are
    joined as ``key=value`` pairs separated by ';'. A key whose value is None is
    skipped, and an empty result yields the '.' placeholder."""
    pairs = []
    for key, value in attributes.items():
        if value is None:
            continue
        pairs.append(f"{_escape_gff_attr(str(key))}={_escape_gff_attr(str(value))}")
    return ";".join(pairs) if pairs else "."


class Feature:
    """A single cross-annotation feature
    `start` and `end` are sorted by smallest to largest and are expected to be 0-based half-open upon input"""

    def __init__(self, fid=None, pid=None, seqid=None, source=None, type=None, start=None, end=None, score=None, strand=None, phase=None,
                 attributes=None, descendants: list = None, parent: "Feature" = None, sequence: str = "", origin_sequence: str = "",
                 ingest: bool = False
                 ):
        # feature ID
        self.fid = fid
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

        if ingest:
            self._ingest()

        if not self.fid:
            raise AttributeError('No ID obtained for feature')


    def _ingest(self):
        """Normalize ``attributes`` and derive ``fid``/``pid`` from them.

        A string ``attributes`` (a raw GFF3 column-9 field) is parsed into a
        dict; ``fid`` is then set from the first present of ID/Id/id and ``pid``
        from the first of Parent/parent. An attribute hit overrides the
        constructor value; a miss leaves it untouched."""
        if isinstance(self.attributes, str):
            self.attributes = _parse_gff_attributes(self.attributes)
        if not self.fid:
            self.fid = _attr_get(self.attributes, ("ID", "Id", "id"))
        if not self.pid:
            self.pid = _attr_get(self.attributes, ("Parent", "parent"))


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
            score, strand, phase, _format_gff_attributes({**id_dict, **self.attributes}),
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
    id_counts = Counter()
    forbidden = set()
    for feature in features:
        fid = _attr_get(feature.attributes, ids)
        if fid in id2feature:
            forbidden.add(fid)
            # make the ID unique by appending a count
            while fid in id2feature:
                id_counts[fid] += 1
                fid = f"{fid}_{id_counts[fid]}"
            for id_ in ids:
                if id_ in feature.attributes:
                    feature.attributes[id_] = fid
            # keep the cached fid in step with the renamed attribute so the
            # FeatureCol index (keyed on fid) sees a distinct feature per
            # multi-segment CDS rather than colliding on the shared id
            feature.fid = fid
        id2feature[fid] = feature

    feature_dict = defaultdict(list)
    for fid, feature in id2feature.items():
        par_id = _attr_get(feature.attributes, parent_ids)
        if par_id:
            # cannot dereplicate
            if par_id in forbidden:
                raise KeyError(f"parent ID {par_id} belongs to multiple features")
            # link features together
            id2feature[par_id].descendants.append(feature)
            feature.parent = id2feature[par_id]
        feature_dict[feature.seqid].append(feature)

    return feature_dict


# GFF `type` tokens grouped into the canonical feature classes FeatureCol
# exposes. RNA is a category: any type ending in "rna"/"transcript" (mRNA,
# tRNA, ncRNA, primary_transcript, ...) resolves to the "rna" class.
_FEATURE_CLASSES = ("gene", "rna", "cds", "exon")


class FeatureCol:
    """An ordered collection of Features linked by Parent <-> ID relationships.

    On construction the supplied Features are hierarchically grouped (see
    `group_features`), populating each Feature's `parent`/`descendants`.
    Features of a canonical class are then reachable both as attributes
    (`.genes`, `.rnas`, `.cds`, `.exons`) and by key (`fl["gene"]`,
    `fl["mRNA"]`); both resolve to the same lists, extracted once by
    `_extract_types`."""

    def __init__(self, features: list = None, group: bool = True):
        self.features = list(features) if features is not None else []
        if group:
            # wire up parent/descendant links in place
            group_features(self.features)
        self._extract_types()

    @staticmethod
    def _class_of(ftype):
        """Return the canonical class ('gene'/'rna'/'cds'/'exon') for a GFF
        `type`, or None if it is not one of them."""
        if not ftype:
            return None
        t = ftype.lower()
        if t.endswith("rna") or t.endswith("transcript"):
            return "rna"
        # gene/cds/exon map to themselves
        if t in _FEATURE_CLASSES:
            return t
        return None

    def _extract_types(self):
        """Bucket features by canonical class and store each list as an
        attribute (`self.genes`, `self.rnas`, `self.cds`, `self.exons`), and
        build the `fid` -> Feature index used by `by_id`/`get`.

        Raises KeyError if two features share an `fid` (after any deduplication
        `group_features` applied), so an ambiguous index is surfaced rather than
        silently collapsed."""
        buckets = {name: [] for name in _FEATURE_CLASSES}
        index = {}
        for feature in self.features:
            cls = self._class_of(feature.type)
            if cls is not None:
                buckets[cls].append(feature)
            if feature.fid in index:
                raise KeyError(f"duplicate feature ID in collection: {feature.fid!r}")
            index[feature.fid] = feature
        self._buckets = buckets
        self._index = index
        self.genes = buckets["gene"]
        self.rnas = buckets["rna"]
        self.cds = buckets["cds"]
        self.exons = buckets["exon"]

    def sort(self):
        """Order features by contig (seqid), then parent-before-descendant,
        then start coordinate; returns self.

        Root features (no parent within this collection) are grouped by seqid
        and ordered by start, then each is emitted immediately before its
        descendants. Every descendant list is sorted by start in place, so the
        depth-first walk here and `Feature.to_gff` share one sibling order. An
        undefined ('.') seqid or start sorts to the end of its group."""
        def _start_key(feature):
            return (feature.start is None, feature.start)

        for feature in self.features:
            feature.descendants.sort(key=_start_key)

        def _walk(feature):
            yield feature
            for descendant in feature.descendants:
                yield from _walk(descendant)

        roots = [f for f in self.features if f.parent is None]
        roots.sort(key=lambda f: (f.seqid is None, f.seqid, _start_key(f)))

        ordered = []
        for root in roots:
            ordered.extend(_walk(root))
        self.features = ordered
        self._extract_types()
        return self

    def roots(self):
        """Return a new FeatureCol of only the root features (those with no
        parent within this collection).

        Each root retains its existing `descendants`/`parent` wiring, so the
        full hierarchy stays reachable through `.descendants`; grouping is
        skipped to preserve those links rather than rebuild them."""
        return FeatureCol([f for f in self.features if f.parent is None], group=False)

    def by_id(self, fid):
        """Return the Feature whose `fid` equals `fid`, raising KeyError if no
        feature in this collection has that ID.

        Only features in `self.features` are indexed; descendants reachable only
        through `.descendants` (e.g. after `roots()`) are not, mirroring how the
        class buckets are populated."""
        return self._index[fid]

    def get(self, fid, default=None):
        """Return the Feature whose `fid` equals `fid`, or `default` if no
        feature in this collection has that ID (non-raising counterpart of
        `by_id`)."""
        return self._index.get(fid, default)

    def _resolve_key(self, key: str):
        """Map an access key to a canonical class: a class name ('gene'), its
        plural ('genes'), or a raw GFF type ('mRNA')."""
        k = key.lower()
        if k in self._buckets:
            return k
        if k.endswith("s") and k[:-1] in self._buckets:
            return k[:-1]
        cls = self._class_of(k)
        if cls is None:
            raise KeyError(f"{key!r} is not a gene/RNA/CDS/exon feature class")
        return cls

    def __getitem__(self, key):
        """Key by feature class ('gene'/'rna'/'cds'/'exon', a plural, or a raw
        GFF type such as 'mRNA') to get that class's list; index or slice into
        the full feature list with an int or slice."""
        if isinstance(key, (int, slice)):
            return self.features[key]
        if isinstance(key, str):
            return self._buckets[self._resolve_key(key)]
        raise TypeError(
            f"FeatureCol keys must be str, int, or slice, not {type(key).__name__}"
        )

    def __iter__(self):
        return iter(self.features)

    def __len__(self):
        return len(self.features)

    def __repr__(self):
        return (f"FeatureCol({len(self.genes)} genes, {len(self.rnas)} RNAs, "
                f"{len(self.cds)} CDS, {len(self.exons)} exons)")