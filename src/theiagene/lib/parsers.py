"""Shared reference-parsing layer for the theiagene commands.
"""

import os
import gzip
import json
import logging
import threading
from collections import defaultdict
from urllib.parse import unquote

import pysam
from Bio import SeqIO

from theiagene.lib.gene_model import Gene, GeneModel, Transcript
from theiagene.lib.query import (
    exact_check,
    substring_check,
    match_query,
    normalize_name,
    sanitize_info_value,
)


logger = logging.getLogger(__name__)


# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
GFF_STRAND = {"+": 1, "-": -1}


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


def iter_gff_features(reference_gff: str, feature_type: str):
    """Yield features of optional ``feature_type`` from a GFF3 file.

    Comment/directive lines, blank lines, malformed (non 9-column) lines and any
    embedded ``##FASTA`` section are skipped.  Each yielded feature is a dict
    with keys ``seqid``, ``type``, ``start`` (0-based), ``end`` (half-open),
    ``strand`` (1/-1/None), ``phase`` and ``attributes`` (a parsed dict).  A
    falsy ``feature_type`` yields every feature (used to read the whole file
    when assembling the parent/child hierarchy)."""
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
            if feature_type:
                if obs_type.lower() != feature_type.lower():
                    continue
            yield {
                "seqid": seqid,
                "source": source,
                "type": obs_type,
                # GFF columns are 1-based, both-inclusive; convert to 0-based, half-open
                "start": int(start) - 1,
                "end": int(end),
                "score": score,
                "strand": GFF_STRAND.get(strand),
                "phase": phase,
                "attributes": parse_gff_attributes(attributes),
            }


def iter_gbff_raw(reference_gbff: str, feature_types={"CDS"}):
    """Yield a :class:`Gene` (coordinates, qualifiers and reference sequence, but
    no resolved identity yet) for every expected feature in a GBFF.

    GBFF features are flat (BioPython does not link parent/child), so each
    yielded gene carries a single default transcript; a coding one is later
    modelled as a ``<gene>_mRNA`` spoof transcript."""
    lower_feature_types = set(x.lower() for x in feature_types)
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_seq = str(record.seq)
            for feature in record.features:
                if feature.type.lower() not in lower_feature_types:
                    continue
                # GenBanks are 1-based; BioPython adjusts to 0-based half-open
                parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
                yield Gene(
                    contig_candidates=[record.id, record.name],
                    strand=feature.location.strand,
                    qualifiers=dict(feature.qualifiers),
                    parts={feature.type.lower(): parts},
                    contig_seq=contig_seq,
                )


def _attr_get(attributes: dict, keys):
    """Return the value of the first present attribute key (accepts case
    variants, e.g. ``ID``/``id`` or ``Parent``/``parent``), or ``None``"""
    for key in keys:
        if key in attributes:
            return attributes[key]
    return None


def _accumulate_qualifiers(qualifiers: dict, attributes: dict) -> None:
    """Merge a feature's attributes into a gene's ``{name: [values]}`` qualifiers,
    de-duplicating values while preserving first-seen order"""
    for name, value in attributes.items():
        values = qualifiers.setdefault(name, [])
        if value not in values:
            values.append(value)


def _feature_id(feature: dict):
    """The feature's ``ID`` attribute (case variants accepted), or ``None``"""
    return _attr_get(feature["attributes"], ("ID", "id", "Id"))


def _resolve_path(feature: dict, id2feature: dict) -> list:
    """Walk a feature's ``Parent`` chain up to its root ancestor.

    Returns ``[feature, ..., root]`` -- every feature reached by following the
    first ``Parent`` link up to the enclosing ``gene`` (or the feature itself when
    it has no resolvable parent).  ``Parent`` may list several ids; the first is
    followed.  A dangling or cyclic link stops the walk at the last resolved
    feature."""
    path = [feature]
    seen = {id(feature)}
    current = feature
    while True:
        parents = _attr_get(current["attributes"], ("Parent", "parent"))
        if not parents:
            break
        parent = id2feature.get(parents.split(",")[0].strip())
        if parent is None or id(parent) in seen:
            break
        seen.add(id(parent))
        path.append(parent)
        current = parent
    return path


def _spoof_transcript_id(root_id):
    """Name for the synthesized transcript wrapping a CDS with no RNA layer"""
    return f"{root_id}_mRNA" if root_id is not None else None


def _transcript_key(path: list, root_id, child_ids: set):
    """The transcript id a feature belongs to, given its resolved ``path``.

    The transcript is the child of the root gene along the path -- the ``mRNA`` in
    a well-formed ``gene -> mRNA -> CDS`` hierarchy.  When a feature hangs directly
    off the gene (or is a bare root) and is not itself an RNA/transcript with
    descendants, its CDS is wrapped in a ``<gene>_mRNA`` spoof transcript so it is
    still translated."""
    feature = path[0]
    if len(path) >= 3:
        # feature -> ... -> transcript -> root: the transcript is the child-of-root
        return _feature_id(path[-2]) or _spoof_transcript_id(root_id)
    # len 2 (feature's parent is the gene) or 1 (feature is itself a root): a
    # feature that is someone's Parent is a real (RNA) transcript; otherwise spoof
    fid = _feature_id(feature)
    if fid is not None and fid in child_ids:
        return fid
    return _spoof_transcript_id(root_id)


def iter_gff_raw(reference_gff: str, fa_dict: dict = None):
    """Yield one :class:`Gene` per gene in a GFF3, partitioned into transcripts.

    Every feature is grouped onto its root ancestor -- the enclosing ``gene`` in a
    well-formed GFF3, resolved by following each feature's ``Parent`` up the chain
    (``gene -> RNA -> CDS``, or ``gene -> CDS`` directly).  Within a gene, each
    feature is filed under its **transcript** (the ``RNA`` child of the gene); a
    CDS whose parent is the gene directly -- or a bare gene-less CDS -- is wrapped
    in a synthesized ``<gene>_mRNA`` **spoof transcript** so it is still modelled.
    Two RNAs under one gene therefore yield two transcripts (two coding models).

    The gene's ``qualifiers`` accumulate the attributes of the root together with
    those of every constituent (so a query matches on ``product`` from a CDS or
    ``gene``/``locus_tag`` from the gene line), each de-duplicated in first-seen
    order; a transcript's own ``qualifiers`` carry only its features' attributes
    (so its product can be resolved per isoform).  A multi-segment CDS sharing an
    ``ID`` under one RNA coalesces into one transcript; a gene with no CDS is
    yielded with empty :attr:`Gene.cds`.

    When ``fa_dict`` (a ``{seqid: SeqRecord}`` mapping) is given, ``contig_seq``
    is filled from it (``None`` when the seqid is absent, so the caller can raise
    with a format-appropriate message); otherwise it is ``None`` (coordinate-only
    use)."""
    # read every feature; index by ID for Parent-chain resolution (gene/RNA ids are
    # unique, multi-segment CDS reuse one id harmlessly -- first occurrence wins)
    features = list(iter_gff_features(reference_gff, None))
    id2feature = {}
    child_ids = set()
    for feature in features:
        fid = _feature_id(feature)
        if fid is not None and fid not in id2feature:
            id2feature[fid] = feature
        parents = _attr_get(feature["attributes"], ("Parent", "parent"))
        if parents:
            for pid in parents.split(","):
                child_ids.add(pid.strip())

    groups = {}
    for feature in features:
        path = _resolve_path(feature, id2feature)
        root = path[-1]
        root_id = _feature_id(root)
        # a root without an ID keys on its object identity so it stays distinct
        group_key = (root["seqid"], root_id if root_id is not None else id(root))
        group = groups.get(group_key)
        if group is None:
            group = {
                "seqid": root["seqid"],
                "strand": root["strand"],
                "root_id": root_id,
                "qualifiers": {},
                "transcripts": {},
                "torder": [],
            }
            # seed with the gene(root) attributes, then layer on each constituent
            _accumulate_qualifiers(group["qualifiers"], root["attributes"])
            groups[group_key] = group
        # gene-level qualifiers union every feature (drives matching)
        _accumulate_qualifiers(group["qualifiers"], feature["attributes"])

        tid = _transcript_key(path, root_id, child_ids)
        transcript = group["transcripts"].get(tid)
        if transcript is None:
            transcript = {"parts": defaultdict(list), "qualifiers": {}, "strand": None}
            group["transcripts"][tid] = transcript
            group["torder"].append(tid)
        transcript["parts"][feature["type"]].append((feature["start"], feature["end"]))
        # per-transcript qualifiers carry only this transcript's features
        _accumulate_qualifiers(transcript["qualifiers"], feature["attributes"])
        if feature["type"].lower() == "cds" and feature["strand"] is not None:
            transcript["strand"] = feature["strand"]

    # materialize each contig once and share it across every gene on it: str(Seq)
    # decodes a fresh copy per call, so deriving it per gene holds as many
    # full-contig copies as there are matched genes (models share by reference)
    contig_seqs = {}

    for group in groups.values():
        seqid = group["seqid"]
        if fa_dict is not None and seqid in fa_dict and seqid not in contig_seqs:
            contig_seqs[seqid] = str(fa_dict[seqid].seq)
        transcripts = {}
        for tid in group["torder"]:
            raw_transcript = group["transcripts"][tid]
            transcripts[tid] = Transcript(
                transcript_id=tid,
                strand=(
                    raw_transcript["strand"]
                    if raw_transcript["strand"] is not None
                    else group["strand"]
                ),
                qualifiers=raw_transcript["qualifiers"],
                parts={k: list(v) for k, v in raw_transcript["parts"].items()},
            )
        yield Gene(
            contig_candidates=[seqid],
            strand=group["strand"],
            qualifiers=group["qualifiers"],
            transcripts=transcripts,
            contig_seq=contig_seqs.get(seqid),
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
    matching is normalization-aware via
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


def _id_qualifiers(feature_qualifier: str) -> list:
    """CDS qualifiers matched for query resolution (query sets may mix product
    names and locus tags), the feature qualifier first"""
    return [feature_qualifier.strip(), "gene", "locus_tag", "protein_id"]


def _register_model(models_by_key: dict, model: GeneModel, key_sources) -> None:
    """Register a model under raw/sanitized/normalized forms of every identifier,
    so it can be recovered from whatever identifier the VCF ``GENE`` field carries"""
    keys = set()
    for ident in key_sources:
        keys.update((ident, sanitize_info_value(ident), normalize_name(ident)))
    for key in keys:
        if not key:
            continue
        if key in models_by_key and models_by_key[key] is not model:
            logger.warning(f"'{key}' recovered multiple times; keeping first")
        else:
            models_by_key[key] = model


def _assemble_model(raw, transcript, matched_query, feature_qualifier, transl_table_override):
    """Generate the GeneModel for one coding transcript of a matched Gene.

    ``product``/``transl_table`` are resolved from the transcript's own qualifiers
    first (so isoforms keep distinct products), falling back to the gene-level
    union.  Only ``raw.contig`` (constant per gene) is stamped on the locus; the
    per-transcript identity is passed to :meth:`GeneModel.from_transcript`.
    Returns (model, gene_id, product)."""
    qualifier = feature_qualifier.strip()
    if transl_table_override is not None:
        transl_table = transl_table_override
    else:
        tt = transcript.qualifiers.get("transl_table") or raw.qualifiers.get("transl_table")
        transl_table = int(tt[0]) if tt else 1
    product_vals = transcript.qualifiers.get(qualifier) or raw.qualifiers.get(qualifier)
    product = product_vals[0] if product_vals else matched_query
    gene_id = normalize_name(matched_query)
    model = GeneModel.from_transcript(
        raw,
        transcript,
        gene_id=gene_id,
        contig=raw.contig,
        product=product,
        transl_table=transl_table,
    )
    return model, gene_id, product


def _register_gene_models(
    models_by_key, raw, matched_query, identifiers, feature_qualifier, transl_table_override
):
    """Build and register a GeneModel for each coding transcript of a matched Gene.

    Each isoform is registered under the gene-level identifiers plus its own
    ``transcript_id``, so distinct isoforms are never lost to keep-first when they
    share a product/gene name."""
    for transcript in raw.coding_transcripts():
        strand = transcript.strand if transcript.strand in (1, -1) else raw.strand
        if strand not in (1, -1):
            logger.warning(
                f"Skipping '{matched_query}' ({transcript.transcript_id}): "
                f"unresolved strand ({strand})"
            )
            continue
        transcript.strand = strand
        model, gene_id, product = _assemble_model(
            raw, transcript, matched_query, feature_qualifier, transl_table_override
        )
        _register_model(
            models_by_key,
            model,
            identifiers + [matched_query, product, gene_id, model.transcript_id],
        )


def build_gene_models_gbff(
    reference_gbff: str,
    contig_names: set,
    query_genes,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GBFF into {lookup_key: GeneModel} for every query gene.

    A CDS is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` qualifiers, so query sets may mix product
    names and locus tags.  A model is registered under many lookup keys (raw,
    sanitized and normalized forms of every identifier) so it can be recovered
    from whatever identifier the VCF ``GENE`` field carries."""
    query_list = list(query_genes)
    id_qualifiers = _id_qualifiers(feature_qualifier)
    models_by_key = {}
    for raw in iter_gbff_raw(reference_gbff):
        matched_query, identifiers = match_identifiers(
            raw.qualifiers, query_list, id_qualifiers,
            exact_match=exact_match, normalize=True,
        )
        if matched_query is None:
            continue
        if not raw.cds:
            logger.warning(f"Skipping '{matched_query}': no CDS coordinates to model")
            continue
        raw.contig = resolve_contig(raw.contig_candidates, contig_names, True, "VCF")
        _register_gene_models(
            models_by_key, raw, matched_query, identifiers,
            feature_qualifier, transl_table_override,
        )
    return models_by_key


def build_gene_models_gff(
    reference_gff: str,
    reference_fa: str,
    contig_names: set,
    query_genes,
    feature_qualifier: str,
    exact_match: bool = False,
    transl_table_override: int = None,
) -> dict:
    """Parse a GFF and genome FA into {lookup_key: GeneModel} for every query gene.

    A gene is matched by the ``feature_qualifier`` value plus its ``gene``,
    ``locus_tag`` and ``protein_id`` attributes, so query sets may mix product
    names and locus tags.  Each of the gene's **coding transcripts** yields its
    own model: multi-exon CDS of one transcript are assembled in translation
    order, while two RNAs under one gene give two models (distinguished by
    ``product``).  Each model is registered under many lookup keys (raw, sanitized
    and normalized forms of every identifier plus its ``transcript_id``) so it can
    be recovered from whatever identifier the VCF ``GENE`` field carries.

    The coding phase column is not applied; as with the GBFF backend, each CDS
    is assumed to begin on a codon boundary."""
    query_list = list(query_genes)
    id_qualifiers = _id_qualifiers(feature_qualifier)
    fa_dict = SeqIO.to_dict(SeqIO.parse(reference_fa, "fasta"))
    models_by_key = {}
    for raw in iter_gff_raw(reference_gff, fa_dict):
        matched_query, identifiers = match_identifiers(
            raw.qualifiers, query_list, id_qualifiers,
            exact_match=exact_match, normalize=True,
        )
        if matched_query is None:
            continue
        if not raw.cds:
            logger.warning(f"Skipping '{matched_query}': no CDS coordinates to model")
            continue
        # check appropriate contig is used for the VCF and available in the FASTA
        contig = resolve_contig(raw.contig_candidates, contig_names, True, "VCF")
        if raw.contig_seq is None:
            raise KeyError(f"{contig} not in reference FASTA")
        raw.contig = contig
        _register_gene_models(
            models_by_key, raw, matched_query, identifiers,
            feature_qualifier, transl_table_override,
        )
    return models_by_key


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
