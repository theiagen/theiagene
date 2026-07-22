"""Shared reference-parsing layer for the theiagene commands.
"""

import logging
from collections import defaultdict
from urllib.parse import unquote

import pysam
from Bio import SeqIO

from theiagene.lib.gene_model import Gene
from theiagene.lib.query import (
    exact_check,
    substring_check,
    match_query,
    sanitize_info_value,
)


logger = logging.getLogger(__name__)


# GFF strand column -> BioPython-style strand integer ('.'/'?' -> None)
GFF_STRAND = {"+": 1, "-": -1}


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


def iter_gbff_raw(reference_gbff: str):
    """Yield a :class:`Gene` (coordinates, qualifiers and reference sequence, but
    no resolved identity yet) for every CDS feature in a GBFF"""
    with open(reference_gbff) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig_seq = str(record.seq)
            for feature in record.features:
                if feature.type.lower() != "cds":
                    continue
                # GenBanks are 1-based; BioPython adjusts to 0-based half-open
                parts = [(int(p.start), int(p.end)) for p in feature.location.parts]
                yield Gene(
                    contig_candidates=[record.id, record.name],
                    strand=feature.location.strand,
                    qualifiers=dict(feature.qualifiers),
                    parts={"CDS": parts},
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


def _resolve_root(feature: dict, id2feature: dict) -> dict:
    """Walk a feature's ``Parent`` chain up to its root ancestor.

    Returns the top-most feature reachable through ``Parent`` links -- the
    enclosing ``gene`` in a well-formed GFF3 (``gene -> RNA -> exon/CDS``), or
    the feature itself when it has no resolvable parent.  ``Parent`` may list
    several ids; the first is followed.  A dangling or cyclic link stops the
    walk at the last resolved feature."""
    current = feature
    seen = {id(current)}
    while True:
        # first parent retrieved
        parents = _attr_get(current["attributes"], ("Parent", "parent"))
        # if no parent, it is the source
        if not parents:
            return current
        parent = id2feature.get(parents.split(",")[0].strip())
        if parent is None or id(parent) in seen:
            return current
        seen.add(id(parent))
        current = parent


def iter_gff_raw(reference_gff: str, fa_dict: dict = None):
    """Yield one :class:`Gene` per gene in a GFF3, assimilating its constituents
    (CDS, exon, RNA, ...) through the parent/child hierarchy.

    Every feature is grouped onto its root ancestor -- the enclosing ``gene`` in a
    well-formed GFF3, resolved by following each feature's ``Parent`` up the chain
    (``gene -> RNA -> CDS``, or ``gene -> CDS`` directly).  The gene is identified
    by the ``ID`` of that root; its coordinate ``parts`` are filed by feature type
    (so :attr:`Gene.cds` yields the CDS segments used for coverage/translation),
    and its ``qualifiers`` accumulate the attributes of the root together with
    those of its constituents (so ``product`` is taken from the CDS while
    ``gene``/``locus_tag`` come from the gene line), each as a de-duplicated list
    in first-seen order.  A feature with no resolvable gene ancestor is grouped on
    its own root ``ID``, so a gene-less GFF still yields one Gene per coding
    feature (multi-segment CDS sharing an ``ID`` still coalesce); a gene with no
    CDS is yielded with an empty :attr:`Gene.cds`.

    When ``fa_dict`` (a ``{seqid: SeqRecord}`` mapping) is given, ``contig_seq``
    is filled from it (``None`` when the seqid is absent, so the caller can raise
    with a format-appropriate message); otherwise it is ``None`` (coordinate-only
    use)."""
    # create a dictionary of GFF features by line
    features = list(iter_gff_features(reference_gff, None))
    # index features by ID for Parent-chain resolution; gene/RNA ids are unique,
    # multi-segment CDS reuse one id harmlessly (first occurrence wins)
    id2feature = {}
    for feature in features:
        fid = _attr_get(feature["attributes"], ("ID", "id", "Id"))
        if fid is not None and fid not in id2feature:
            id2feature[fid] = feature

    groups = {}
    for feature in features:
        root = _resolve_root(feature, id2feature)
        root_id = _attr_get(root["attributes"], ("ID", "id"))
        # a root without an ID keys on its object identity so it stays distinct
        group_key = (root["seqid"], root_id if root_id is not None else id(root))
        group = groups.get(group_key)
        if group is None:
            group = {
                "seqid": root["seqid"],
                "strand": root["strand"],
                "qualifiers": {},
                "parts": defaultdict(list),
            }
            # seed with the gene(root) attributes, then layer on each constituent
            _accumulate_qualifiers(group["qualifiers"], root["attributes"])
            groups[group_key] = group
        group["parts"][feature["type"]].append((feature["start"], feature["end"]))
        _accumulate_qualifiers(group["qualifiers"], feature["attributes"])

    for group in groups.values():
        seqid = group["seqid"]
        yield Gene(
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


def flatten_coords_by_contig(genes, full_range: bool = False) -> dict:
    """Flatten Gene objects to {<CONTIG>: [(START, END, GENE_ID), ...]} for interval
    overlap testing.  If full_range is True, each gene is collapsed to a single
    (genomic_start, genomic_end) range spanning all of its parts; otherwise every
    part is emitted separately"""
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


def extract_vcf_genes(vcffile: str, genes, output_vcf: str) -> int:
    """Filter a VCF to variants overlapping query gene coordinates, annotating the
    overlapping gene name(s) in a GENE INFO field. Returns the count of written records"""
    vcf_in = pysam.VariantFile(vcffile)
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
