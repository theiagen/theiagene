"""Gene locus model: :class:`Transcript`, :class:`Gene` and :class:`GeneModel`.

A :class:`Transcript` is one coding unit -- its genomic coordinate segments
(0-based, half-open), filed by feature type and ordered into translation
(5'->3') order by ``strand``.  A :class:`Gene` is the *locus*: metadata plus a
collection of transcripts, exposing coverage-facing **union** views over them.
A :class:`GeneModel` is a :class:`Transcript` that additionally carries its
reference sequence and, from it, the coding sequence and four derived sequences;
one is built per coding transcript of a gene.
"""

import logging

from theiagene.lib.sequence import complement, reverse_complement, translate


logger = logging.getLogger(__name__)


# sentinel key for the anonymous default transcript a Gene creates when parts are
# added without a transcript id (a bare ``parts=`` kwarg, ``add_part`` with no
# ``transcript_id``, BED rows, or a GBFF CDS)
_DEFAULT_TRANSCRIPT = object()


class Transcript:
    """Coordinates for a single coding unit (one transcript / one CDS chain).

    ``parts`` is a ``{feature_type: [(start, end), ...]}`` dict -- CDS segments
    live under a ``"CDS"`` key (see :attr:`cds`); the derived properties order
    those into translation (5'->3') order using ``strand`` (1/-1/None).
    ``qualifiers`` holds this transcript's own parsed feature attributes
    (``{name: [values]}``), used to resolve its product."""

    def __init__(
        self,
        transcript_id=None,
        strand=None,
        transl_table=1,
        product=None,
        parts=None,
        qualifiers=None,
    ):
        self.transcript_id = transcript_id
        self.strand = strand
        self.transl_table = transl_table
        self.product = product
        # parts must be a {feature_type: [(start, end), ...]} dict (0-based,
        # half-open); CDS segments live under a "CDS" key (see the cds property)
        self.parts = {}
        if parts:
            if not isinstance(parts, dict):
                raise TypeError(
                    "'parts' must be a {feature_type: [(start, end), ...]} "
                    f"dict, not {type(parts).__name__}"
                )
            for feature, segments in parts.items():
                self.parts[feature] = [(int(s), int(e)) for s, e in segments]
        self.qualifiers = qualifiers if qualifiers is not None else {}

    def segments(self, feature_type: str) -> list:
        """Coordinate segments filed under ``feature_type`` (case-insensitive key
        match), or ``[]`` when the transcript carries no part of that type
        (e.g. ``"CDS"``, ``"exon"``); :attr:`cds` is the ``"CDS"`` specialisation."""
        target = feature_type.lower()
        for key, segments in self.parts.items():
            if key.lower() == target:
                return segments
        return []

    @property
    def cds(self) -> list:
        """CDS coordinate segments (case-insensitive key lookup); ``[]`` if none"""
        return self.segments("CDS")

    def add_part(self, start, end, feature: str = "CDS") -> None:
        """Append a genomic ``(start, end)`` segment (0-based, half-open) under
        ``feature`` (the CDS coordinate list by default)"""
        self.parts.setdefault(feature, []).append(
            tuple(sorted((int(start), int(end))))
        )

    @property
    def genomic_start(self):
        """Left-most genomic coordinate across the CDS parts (None if empty)"""
        return min((s for s, _ in self.cds), default=None)

    @property
    def genomic_end(self):
        """Right-most genomic coordinate across the CDS parts (None if empty)"""
        return max((e for _, e in self.cds), default=None)

    @property
    def cds_positions(self):
        """Every coding-base genomic position in translation (5'->3') order"""
        positions = []
        for start, end in sorted(self.cds):
            positions.extend(range(start, end))
        if self.strand == -1:
            positions.reverse()
        return positions


class Gene:
    """A gene locus: metadata plus its transcripts.

    Coordinates are 0-based, half-open, on the reference contig.  A gene holds a
    ``{transcript_id: Transcript}`` mapping (see :attr:`transcripts`); a
    :class:`GeneModel` is built per coding transcript.  The **union views**
    (:attr:`parts`, :meth:`segments`, :attr:`cds`, :attr:`genomic_start`/
    :attr:`genomic_end`, :attr:`cds_positions`) aggregate across all transcripts
    and back the coverage/overlap paths and the legacy single-coding-unit
    construction; per-isoform translation uses each transcript's own coordinates.
    :attr:`cds_positions` is only meaningful for a single-transcript gene.

    ``qualifiers`` holds the parsed feature attributes (``{name: [values]}``) used
    to match the gene against a query set, ``contig_candidates`` the seqid/name(s)
    a resolved ``contig`` is chosen from, and ``contig_seq`` the reference
    sequence when available (required to build a :class:`GeneModel`)."""

    def __init__(
        self,
        gene_id=None,
        contig=None,
        strand=None,
        transl_table=1,
        product=None,
        parts=None,
        qualifiers=None,
        contig_seq=None,
        contig_candidates=None,  # check multiple contig names for cross-format compatibility
        transcripts=None,
    ):
        self.gene_id = gene_id
        self.contig = contig
        self.strand = strand
        self.transl_table = transl_table
        self.product = product if product is not None else gene_id
        self.qualifiers = qualifiers if qualifiers is not None else {}
        self.contig_seq = contig_seq
        self.contig_candidates = (
            list(contig_candidates)
            if contig_candidates is not None
            else ([contig] if contig is not None else [])
        )
        self.transcripts = {}
        if transcripts:
            for tid, transcript in transcripts.items():
                self.transcripts[tid] = transcript
        if parts:
            # a bare parts dict seeds the anonymous default transcript, so the
            # legacy single-coding-unit Gene(parts=...) construction keeps working
            if not isinstance(parts, dict):
                raise TypeError(
                    "Gene 'parts' must be a {feature_type: [(start, end), ...]} "
                    f"dict, not {type(parts).__name__}"
                )
            default = self._default_transcript()
            for feature, segments in parts.items():
                default.parts.setdefault(feature, []).extend(
                    (int(s), int(e)) for s, e in segments
                )

    def _default_transcript(self) -> Transcript:
        """Get (creating if needed) the anonymous default transcript"""
        transcript = self.transcripts.get(_DEFAULT_TRANSCRIPT)
        if transcript is None:
            transcript = Transcript(
                transcript_id=self.gene_id,
                strand=self.strand,
                transl_table=self.transl_table,
                product=self.product,
            )
            self.transcripts[_DEFAULT_TRANSCRIPT] = transcript
        return transcript

    def add_transcript(self, transcript: Transcript) -> None:
        """Attach a fully-built :class:`Transcript` under its ``transcript_id``"""
        self.transcripts[transcript.transcript_id] = transcript

    def add_part(self, start, end, feature: str = "CDS", transcript_id=None) -> None:
        """Append a genomic ``(start, end)`` segment (0-based, half-open) under
        ``feature`` to a transcript (the default transcript when ``transcript_id``
        is None, else the named one, created on first use)"""
        if transcript_id is None:
            transcript = self._default_transcript()
        else:
            transcript = self.transcripts.get(transcript_id)
            if transcript is None:
                transcript = Transcript(
                    transcript_id=transcript_id,
                    strand=self.strand,
                    transl_table=self.transl_table,
                    product=self.product,
                )
                self.transcripts[transcript_id] = transcript
        transcript.add_part(start, end, feature)

    def coding_transcripts(self) -> list:
        """The transcripts carrying CDS coordinates, in insertion order"""
        return [t for t in self.transcripts.values() if t.cds]

    @property
    def parts(self) -> dict:
        """Union of every transcript's parts as a ``{feature_type: [(start, end), ...]}`` dict"""
        merged = {}
        for transcript in self.transcripts.values():
            for feature, segments in transcript.parts.items():
                merged.setdefault(feature, []).extend(segments)
        return merged

    def segments(self, feature_type: str) -> list:
        """Union of the ``feature_type`` segments across all transcripts
        (case-insensitive key match); ``[]`` when none carry that type"""
        target = feature_type.lower()
        collected = []
        for transcript in self.transcripts.values():
            for key, segments in transcript.parts.items():
                if key.lower() == target:
                    collected.extend(segments)
        return collected

    @property
    def cds(self) -> list:
        """Union of the CDS segments across all transcripts; ``[]`` if none"""
        return self.segments("CDS")

    @property
    def genomic_start(self):
        """Left-most genomic coordinate across the union CDS parts (None if empty)"""
        return min((s for s, _ in self.cds), default=None)

    @property
    def genomic_end(self):
        """Right-most genomic coordinate across the union CDS parts (None if empty)"""
        return max((e for _, e in self.cds), default=None)

    @property
    def cds_positions(self):
        """Every coding-base genomic position (5'->3') across the union CDS.

        Only meaningful for a single-transcript gene; per-isoform translation
        uses each :class:`GeneModel`'s own :attr:`Transcript.cds_positions`."""
        positions = []
        for start, end in sorted(self.cds):
            positions.extend(range(start, end))
        if self.strand == -1:
            positions.reverse()
        return positions


class GeneModel(Transcript):
    """A :class:`Transcript` plus its reference sequence and derived sequences.

    :meth:`finalize` translates the transcript's coordinates against a contig
    sequence into the coding sequence (``ref_coding``, in 5'->3' translation
    order, reverse-complemented for minus-strand) and its ``pos2cds`` index --
    both used by the variant placement/codon math in ``theiagene.lib.variant`` --
    together with four sequence attributes:

    ``protein``      translation of the spliced coding sequence
    ``rna``          the spliced coding sequence as RNA (coding strand, T->U)
    ``dna``          the full transcript span, introns included, on the coding strand
    ``revcomp_dna``  the reverse complement of ``dna`` (template strand)

    Beyond the transcript's coordinates it carries the gene-level identity the
    downstream annotation needs (``gene_id``, ``contig``, ``contig_seq``,
    ``contig_candidates``)."""

    def __init__(
        self,
        gene_id=None,
        contig=None,
        strand=None,
        transl_table=1,
        product=None,
        parts=None,
        qualifiers=None,
        contig_seq=None,
        contig_candidates=None,
        transcript_id=None,
    ):
        super().__init__(
            transcript_id=transcript_id if transcript_id is not None else gene_id,
            strand=strand,
            transl_table=transl_table,
            product=product if product is not None else gene_id,
            parts=parts,
            qualifiers=qualifiers,
        )
        self.gene_id = gene_id
        self.contig = contig
        self.contig_seq = contig_seq
        self.contig_candidates = (
            list(contig_candidates)
            if contig_candidates is not None
            else ([contig] if contig is not None else [])
        )
        self.ref_coding = ""
        self.pos2cds = {}
        self.protein = ""
        self.rna = ""
        self.dna = ""
        self.revcomp_dna = ""

    @classmethod
    def from_gene(cls, gene: Gene) -> "GeneModel":
        """Generate a finalized :class:`GeneModel` from a :class:`Gene`.

        Builds from the gene's **union** CDS coordinates (its :attr:`Gene.parts`),
        requiring CDS coordinates and a reference :attr:`~Gene.contig_seq`, raising
        ``ValueError`` if either is missing.  For a multi-transcript gene this
        models the *union* of all CDS (legacy single-coding-unit semantics) and
        warns; per-isoform models come from :meth:`from_transcript`."""
        if not gene.cds:
            raise ValueError(f"cannot model '{gene.gene_id}': no CDS coordinates")
        if gene.contig_seq is None:
            raise ValueError(f"cannot model '{gene.gene_id}': no reference sequence")
        coding = gene.coding_transcripts()
        if len(coding) > 1:
            logger.warning(
                f"'{gene.gene_id}': modeling the union of {len(coding)} coding "
                "transcripts; use from_transcript for per-isoform models"
            )
        model = cls(
            gene_id=gene.gene_id,
            contig=gene.contig,
            strand=gene.strand,
            transl_table=gene.transl_table,
            product=gene.product,
            parts=gene.parts,
            qualifiers=gene.qualifiers,
            contig_seq=gene.contig_seq,
            contig_candidates=gene.contig_candidates,
        )
        model.finalize()
        return model

    @classmethod
    def from_transcript(
        cls,
        gene: Gene,
        transcript: Transcript,
        *,
        gene_id=None,
        contig=None,
        product=None,
        transl_table=None,
        strand=None,
    ) -> "GeneModel":
        """Generate a finalized :class:`GeneModel` for one transcript of a gene.

        Identity fields are taken from the explicit keyword when given, else the
        transcript, else the gene.  A transcript with no id is named
        ``<gene_id>_mRNA``.  Requires the transcript to carry CDS coordinates and
        the gene to carry a reference sequence."""
        if not transcript.cds:
            raise ValueError(
                f"cannot model transcript '{transcript.transcript_id}': no CDS coordinates"
            )
        if gene.contig_seq is None:
            raise ValueError(f"cannot model '{gene.gene_id}': no reference sequence")
        gid = gene_id if gene_id is not None else gene.gene_id
        tid = transcript.transcript_id
        if tid is None:
            tid = f"{gid}_mRNA" if gid else None
        resolved_strand = strand
        if resolved_strand is None:
            resolved_strand = (
                transcript.strand if transcript.strand in (1, -1) else gene.strand
            )
        model = cls(
            gene_id=gid,
            contig=contig if contig is not None else gene.contig,
            strand=resolved_strand,
            transl_table=(
                transl_table
                if transl_table is not None
                else (transcript.transl_table or gene.transl_table)
            ),
            product=(
                product if product is not None else (transcript.product or gene.product)
            ),
            parts=transcript.parts,
            qualifiers=(transcript.qualifiers or gene.qualifiers),
            contig_seq=gene.contig_seq,
            contig_candidates=gene.contig_candidates,
            transcript_id=tid,
        )
        model.finalize()
        return model

    def finalize(self, contig_seq: str = None) -> None:
        """Build the coding sequence, position index and derived sequences.

        Uses ``contig_seq`` when given, otherwise the stored
        :attr:`~GeneModel.contig_seq`; raises ``ValueError`` when neither is set."""
        if contig_seq is None:
            contig_seq = self.contig_seq
        if contig_seq is None:
            raise ValueError(f"cannot finalize '{self.gene_id}': no reference sequence")
        self.contig_seq = contig_seq
        positions = self.cds_positions
        coding = "".join(contig_seq[p].upper() for p in positions)
        if self.strand == -1:
            # positions are already reversed; complement each base to get revcomp
            coding = complement(coding)
        self.ref_coding = coding
        self.pos2cds = {g: i for i, g in enumerate(positions)}
        self.protein = translate(self.ref_coding, self.transl_table)
        self.rna = coding.replace("T", "U")
        # full transcript span including introns, on the coding strand
        span = contig_seq[self.genomic_start : self.genomic_end].upper()
        if self.strand == -1:
            span = reverse_complement(span)
        self.dna = span
        self.revcomp_dna = reverse_complement(span)

    def codon(self, codon_number: int) -> str:
        """Return the reference codon (1-based codon number) or '' if incomplete"""
        start = (codon_number - 1) * 3
        codon = self.ref_coding[start : start + 3]
        return codon if len(codon) == 3 else ""

    def aa_at(self, codon_number: int) -> str:
        """Return the reference amino acid (one letter) at a 1-based codon"""
        codon = self.codon(codon_number)
        return translate(codon, self.transl_table) if codon else "X"
