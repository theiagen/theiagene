"""Shared reference-parsing layer for the theiagene commands.
"""

import os
import gzip
import json
import logging
import threading
from collections import Counter
from urllib.parse import unquote

import pysam

from theiagene.lib.feature import Feature, FeatureCol


logger = logging.getLogger(__name__)


# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
GFF_STRAND = {"+": 1, "-": -1}
# strand integer -> GFF strand column (anything else, e.g. None, serializes to '.')
GFF_STRAND_REVERSE = {1: "+", -1: "-"}


def write_json(filename: str, data: dict) -> None:
    """Write a JSON file compatible with WDL"""
    with open(filename, "w") as f:
        if data:
            json.dump(data, f, indent=4)
        else:
            # spoof Cromwell (Terra WDL)
            f.write('{"": 0}')


def parse_gff_attributes(attributes: str, field_delimiter: str = ";", value_delimiter: str = "=") -> dict:
    """Parse a GFF3 attribute string into a {key: value} dict.

    An empty attribute column yields an empty dict; a malformed field 
    (one lacking the value delimiter) raises ``ValueError``.  
    Duplicate keys keep the last occurrence."""
    stripped = attributes.strip().strip(";")
    if not stripped:
        return {}
    parsed = {}
    for field in stripped.split(field_delimiter):
        field = field.strip()
        try:
            key, value = field.split(value_delimiter, 1)
        except ValueError:
            raise ValueError(f"unexpected attributes field: {field}")
        if key.strip() in parsed:
            logger.warning(f"Duplicate key in GFF3 attributes: {key.strip()}")
        parsed[key.strip()] = unquote(value.strip())
    return parsed


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


def format_gff_attributes(attributes: dict, field_delimiter: str = ";", value_delimiter: str = "=") -> str:
    """Serialize a parsed ``{key: value}`` attribute dict back to a GFF3 column-9
    string (the inverse of :func:`parse_gff_attributes`).

    Reserved characters in keys and values are percent-encoded; an empty dict
    yields ``.`` (the GFF3 "no attributes" placeholder)."""
    if not attributes:
        return "."
    return field_delimiter.join(
        f"{_gff_escape(str(key))}{value_delimiter}{_gff_escape(str(value))}"
        for key, value in attributes.items()
    )


def iter_gff_features(reference_gff: str):
    """Yield a Feature class from a GFF3 file.

    A ``.gz`` suffix routes the file through ``gzip`` (text mode), so a compressed
    GFF3 enters the same parsing loop as a plain-text one."""
    opener = gzip.open if reference_gff.endswith(".gz") else open
    # a per-type running count backs the IDs synthesized for records that carry
    # no ID of their own (spec-legal for childless, single-line CDS/exon rows)
    type_counts = Counter()
    with opener(reference_gff, "rt") as handle:
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
            seqid, source, obs_type, start, end, score, strand, phase, raw_attributes = fields
            feature = Feature(
                seqid=seqid,
                source=source,
                type=obs_type,
                # GFF columns are 1-based, both-inclusive; convert to 0-based, half-open
                start=int(start) - 1,
                end=int(end),
                score=score,
                strand=strand,
                phase=phase,
                # let the Feature parse the raw column-9 string and derive fid/pid
                attributes=raw_attributes,
                ingest=True,
            )
            if not feature.fid:
                # synthesize a stable ID so a spec-legal ID-less record still parses
                type_key = (feature.type or "").lower()
                type_counts[type_key] += 1
                feature.synthesize_id(type_counts[type_key])
            yield feature


def assimilate_gff(gff: str) -> FeatureCol:
    """Group every GFF3 record onto its root ancestor Feature via ``Parent``/``ID`` links.

    Returns a :class:`FeatureCol`, which wires up the parent/descendant hierarchy
    on construction. Records sharing an id -- the multi-segment CDS of one RNA,
    which GFF3 permits to repeat a single ``ID`` -- are de-replicated by
    suffixing a count (``cds-1``, ``cds-1_1``, ``cds-1_2``, ...), so each segment
    is indexed distinctly rather than collapsing onto the first occurrence. A
    ``Parent`` pointing at such a repeated id is ambiguous and raises
    ``KeyError``."""
    return FeatureCol(iter_gff_features(gff))


def import_bam(bamfile: str) -> pysam.AlignmentFile:
    """Open a BAM (indexing it first if needed)."""
    imported_bam = pysam.AlignmentFile(bamfile)
    # generate an index if it does not exist
    if not imported_bam.has_index():
        logger.debug("Generating BAM index")
        pysam.index(bamfile)
        imported_bam = pysam.AlignmentFile(bamfile)

    return imported_bam


def iter_bed_rows(bedfile: str, min_columns: int = 4):
    """Yield the whitespace-split columns of each BED data row.

    Comment lines (starting with ``#``) and blank/whitespace-only lines are
    skipped; a data row carrying fewer than ``min_columns`` columns raises
    ``ValueError`` naming the file and 1-based line number. This is the single
    place the BED readers agree on what a data row is, so comment/blank handling
    and the column-count guard live here rather than being re-derived (or
    forgotten) at each call site."""
    with open(bedfile) as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith("#") or not line.strip():
                continue
            data = line.split()
            if len(data) < min_columns:
                raise ValueError(
                    f"{bedfile}:{lineno}: BED row has {len(data)} columns, "
                    f"need >= {min_columns}"
                )
            yield data


def _is_bgzf(path: str) -> bool:
    """True when ``path`` starts with the gzip/BGZF magic bytes.

    A plain-text VCF starts with ``##`` (0x23 0x23), so this cleanly separates a
    compressed input (which pysam probes for an index) from an uncompressed one
    (which neither supports nor needs an index)."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _clean_gq_scientific(line: str) -> str:
    """Rewrite GQ subfields written in scientific/float notation as plain
    integers on a single VCF data line.

    htslib parses GQ per its Integer header type and rejects values like
    ``3.5e+01`` or ``35.0`` — pysam raises even while merely iterating records,
    so the fix has to happen in the text before pysam ever sees the value.
    Coercion is a plain ``int(float(x))``; if a value can't be coerced the
    error is left to propagate, which is fine because pysam couldn't consume
    that value either. Missing values (``.``) are passed through untouched."""
    cols = line.rstrip("\n").split("\t")
    # need FORMAT (col 9) plus at least one sample column (col 10+)
    if len(cols) < 10:
        return line
    fmt = cols[8].split(":")
    if "GQ" not in fmt:
        return line
    gq = fmt.index("GQ")
    for i in range(9, len(cols)):
        sub = cols[i].split(":")
        if gq < len(sub) and sub[gq] != ".":
            sub[gq] = str(int(float(sub[gq])))
            cols[i] = ":".join(sub)
    return "\t".join(cols) + "\n"


def import_vcf(vcffile: str) -> pysam.VariantFile:
    """Open a VCF/BCF, first cleaning GQ values written in scientific/float
    notation into plain integers.

    Some callers emit GQ as e.g. ``3.5e+01``; pysam parses GQ per its Integer
    header type and raises on such values — even iterating the records fails —
    so the file cannot be scrubbed through pysam itself. Instead the VCF text is
    read (decompressing BGZF), rewritten line by line, and streamed straight
    into ``pysam.VariantFile`` through an OS pipe, avoiding an intermediate file
    on disk. pysam needs a real file descriptor (a ``BytesIO`` has no
    ``fileno``), hence the pipe. Cleaning happens eagerly here, so an
    uncoercible GQ raises from this call rather than mid-iteration."""
    opener = gzip.open if _is_bgzf(vcffile) else open
    with opener(vcffile, "rt") as src:
        cleaned = "".join(
            line if line.startswith("#") else _clean_gq_scientific(line)
            for line in src
        ).encode()

    read_fd, write_fd = os.pipe()

    def _feed():
        # the pipe buffer is small, so pump the cleaned bytes from a thread
        # while pysam drains the read end; only raw I/O happens here (all
        # coercion already ran above), so this thread cannot raise a value error
        with os.fdopen(write_fd, "wb") as dst:
            try:
                dst.write(cleaned)
            except BrokenPipeError:
                pass  # reader closed early; nothing left to do

    threading.Thread(target=_feed, daemon=True).start()
    return pysam.VariantFile(os.fdopen(read_fd, "rb"))