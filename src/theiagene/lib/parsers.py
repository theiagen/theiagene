"""Shared reference-parsing layer for the theiagene commands.
"""

import os
import gzip
import json
import logging
import threading
from collections import defaultdict, Counter
from urllib.parse import unquote

import pysam
from Bio import SeqIO

from theiagene.lib.feature import Feature
from theiagene.lib.query import (
    exact_check,
    substring_check,
    sanitize_info_value,
)


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


def iter_gff_features(reference_gff: str, id_qualifiers: list = ["id", "ID", "Id"], parent_qualifiers: list = ["parent", "Parent"]):
    """Yield a Feature class from a GFF3 file."""
    with open(reference_gff) as handle:
        for line in handle:
            fid = None
            pid = None
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
            attributes = parse_gff_attributes(raw_attributes)
            for key in id_qualifiers:
                if key in attributes:
                    fid = attributes[key]
                    break
            if not fid:
                raise ValueError(f"none of {id_qualifiers} found in record attributes") 
            del attributes[key]
            for key in parent_qualifiers:
                if key in attributes:
                    pid = attributes[key]
                    break
            if pid:
                del attributes[key]
            yield Feature(
                fid=fid,
                pid=pid,
                seqid=seqid,
                source=source,
                type=obs_type,
                # GFF columns are 1-based, both-inclusive; convert to 0-based, half-open
                start=int(start) - 1,
                end=int(end),
                score=score,
                strand=strand,
                phase=phase,
                attributes=parse_gff_attributes(attributes),
            )


def _assimilate_gbff_features(features: list) -> None:
    """Populate the ``parent``/``descendants`` links of one record's features in
    place, inferring the hierarchy GenBank leaves implicit.

    GenBank features carry no explicit ``Parent``/``ID`` links (unlike GFF3), so
    nesting is derived two ways.  A CDS/exon is linked to its parent ``*RNA`` by
    a shared ``transcript_id`` qualifier whenever both carry one -- this holds
    even when the file order interleaves isoforms.  Failing that, parse order
    decides: a ``gene`` opens a new locus, an ``*RNA`` nests under the most
    recent gene, and a CDS/exon (or any other subfeature) nests under the most
    recent RNA, falling back to the most recent gene.  Features preceding the
    first gene (e.g. the record-spanning ``source``) stay parentless roots."""
    # map each RNA transcript to its Feature so a CDS/exon carrying the same
    # transcript_id links straight to the right isoform regardless of file order
    transcript2rna = {
        feature.attributes["transcript_id"]: feature
        for feature in features
        if feature.type.endswith("RNA") and "transcript_id" in feature.attributes
    }

    super_gene = False
    current_gene = None
    current_rna = None
    for feature in features:
        if feature.type == "gene":
            current_gene, current_rna = feature, None
            continue
        # subsupergene feature
        if super_gene:
            parent = super_gene
            super_gene = False
        # supergene feature (reset if we're still super)
        if not current_gene:
            super_gene = feature
        if feature.type.endswith("RNA"):
            parent = current_gene
            current_rna = feature
        else:
            # CDS/exon (or any other subfeature): prefer the transcript_id link,
            # else nest under the most recent RNA, else the most recent gene
            transcript_id = feature.attributes.get("transcript_id")
            parent = transcript2rna.get(transcript_id) or current_rna or current_gene
        if parent is not None:
            feature.parent = parent
            feature.pid = parent.id
            parent.descendants.append(feature)


def assimilate_gbff_raw(gbff: str, use_id: bool = True, feature_qualifier: str = "locus_tag") -> list:
    """Parse a GenBank flat file into a flat list of hierarchically assimilated
    :class:`Feature` objects.

    Every feature across every record is materialized as a :class:`Feature`.  The
    seqid is taken from ``record.id`` (an accession.version) when ``use_id`` is
    True, else from ``record.name`` (the LOCUS name).  BioPython feature
    locations are already 0-based, half-open, matching Feature's internal
    convention, so no coordinate conversion is applied.  Every feature qualifier
    present is carried into ``attributes``; BioPython stores each qualifier value
    as a list, so a lone value is unwrapped and multiple values are joined,
    keeping the scalar ``{key: value}`` convention Feature shares with the GFF
    parser.  The full record sequence is passed as ``origin_sequence`` so each
    Feature slices out (and stores) only its own forward-strand span.

    Features are then linked into parent/child trees *within each record* (see
    :func:`_assimilate_gbff_features`): every returned Feature has its ``parent``
    and ``descendants`` populated.  The returned list holds every feature in file
    order, so a root feature precedes its descendants -- mirroring the flat,
    already-linked list :func:`group_features` produces for GFF3."""
    features = []
    fids = Counter()
    with open(gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seqid = record.id if use_id else record.name
            origin_sequence = str(record.seq)
            record_features = []
            for feature in record.features:
                attributes = {
                    key: value[0] if len(value) == 1 else ", ".join(map(str, value))
                    for key, value in feature.qualifiers.items()
                }
                fid = attributes[feature_qualifier]
                # add 1 to the counter (1-indexed IDs)
                fids[fid] += 1
                record_features.append(Feature(
                    fid=attributes[feature_qualifier] + f"_{feature.type}{fids[fid]}",
                    seqid=seqid,
                    type=feature.type,
                    start=int(feature.location.start),
                    end=int(feature.location.end),
                    strand=feature.location.strand,
                    attributes=attributes,
                    origin_sequence=origin_sequence,
                ))
            _assimilate_gbff_features(record_features)
            features.extend(record_features)
        
    return features


def assimilate_gbff(gbff: str) -> list:
    """Group every GBFF record onto its root ancestor Feature via ``Parent``/``ID`` links."""
    # read every feature; index by ID for Parent-chain resolution (gene/RNA ids are
    # unique, multi-segment CDS reuse one id harmlessly -- first occurrence wins)
    features = assimilate_gbff_raw(gbff)
    features_dict = defaultdict(list)
    for feature in features:
        features_dict[feature.seqid].append(feature)
    return features_dict


def assimilate_gff(gff: str) -> list:
    """Group every GFF3 record onto its root ancestor Feature via ``Parent``/``ID`` links."""
    # read every feature; index by ID for Parent-chain resolution (gene/RNA ids are
    # unique, multi-segment CDS reuse one id harmlessly -- first occurrence wins)
    features = list(iter_gff_features(gff, None))
    features_dict = group_features(features)
    return features_dict


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


def import_bam(bamfile: str, ambiguous_contig: bool) -> pysam.AlignmentFile:
    """Open a BAM (indexing it first if needed).

    Raises ``ValueError`` when ``ambiguous_contig`` is requested but the
    reference has more than one contig (the single-contig assumption fails)."""
    imported_bam = pysam.AlignmentFile(bamfile)
    # generate an index if it does not exist
    if not imported_bam.has_index():
        logger.debug("Generating BAM index")
        pysam.index(bamfile)
        imported_bam = pysam.AlignmentFile(bamfile)

    # determine if import is compatible with a single contig reference
    contig_names = imported_bam.references
    if ambiguous_contig:
        # can't apply ambiguous contig approach if there are multiple contigs
        if len(contig_names) > 1:
            raise ValueError(
                "can't use ambiguous_contig coordinates when there are multiple contigs in the reference"
            )
    return imported_bam


def parse_bed_genes(
    bedfile: str,
    query_list,
    exact_match: bool,
    contig_names,
    require: bool = True,
    source_label: str = "BAM",
    feature_type: str = "CDS",
) -> list:
    """Parse a BED file into :class:`Gene` objects (coverage-only format).

    Rows sharing a name on the same contig accumulate as multiple parts, matching
    the multi-segment handling of the GBFF/GFF parsers.  A BED file carries no
    feature type of its own, so its regions are filed under ``feature_type`` (the
    type the coverage caller quantifies) so they are read back consistently."""
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
            gene.add_part(int(data[1]), int(data[2]), feature=feature_type)
    return [genes[key] for key in order]


def flatten_coords_by_contig(genes, full_range: bool = False) -> dict:
    """Flatten Gene objects to {<CONTIG>: [(START, END, GENE_ID), ...]} for interval
    overlap testing.  If full_range is True, each gene is collapsed to a single
    (genomic_start, genomic_end) range spanning all of its parts; otherwise every
    CDS part is emitted separately"""
    contig2ranges = defaultdict(list)
    for gene in genes:
        if full_range:
            # collapse all parts of a gene into a single spanning range
            contig2ranges[gene.contig].append(
                (gene.genomic_start, gene.genomic_end, gene.gene_id)
            )
        else:
            for start, end in gene.cds:
                contig2ranges[gene.contig].append((start, end, gene.gene_id))
    return contig2ranges


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


def extract_vcf_genes(vcf_in, genes, output_vcf: str) -> int:
    """Filter a VCF to variants overlapping query gene coordinates, annotating the
    overlapping gene name(s) in a GENE INFO field. Returns the count of written records"""
    # define the INFO field used to annotate the overlapping gene name(s)
    if "GENE" not in vcf_in.header.info:
        vcf_in.header.info.add(
            "GENE",
            ".",
            "String",
            "Query gene(s) whose extracted coordinate range overlaps this variant",
        )
    vcf_out = pysam.VariantFile(output_vcf, "w", header=vcf_in.header)

    # {<CONTIG>: [(START, END, GENE_ID), ...]} (0-based, half-open coordinates)
    contig2ranges = flatten_coords_by_contig(genes)

    written = 0
    for record in vcf_in:
        ranges = contig2ranges.get(record.contig)
        if not ranges:
            continue
        # pysam VariantRecord coordinates are 0-based, half-open (record.start, record.stop)
        genes = set(
            gene
            for start, end, gene in ranges
            if record.start < end and record.stop > start
        )
        if genes:
            # dedupe while preserving order and sanitize for the INFO field
            clean_genes = sorted([
                sanitize_info_value(gene) for gene in dict.fromkeys(genes)
            ])
            record.info["GENE"] = clean_genes
            vcf_out.write(record)
            written += 1

    vcf_out.close()
    vcf_in.close()
    logger.debug(f"Wrote {written} gene-overlapping variant(s) to {output_vcf}")
    return written