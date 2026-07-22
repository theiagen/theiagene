"""Gene coordinate model (:class:`Gene`) and its sequence-bearing subclass.

A :class:`Gene` is the shared coordinate currency of theiagene: it captures the
metadata and genomic segments of one query gene, and optionally the reference
sequence it sits on.  The coordinate parsers in ``theiagene.lib.parsers`` produce
:class:`Gene` objects directly -- ``gene_coverage`` uses them as-is (it only
needs coordinates), while ``variant_annotation`` upgrades each into the
:class:`GeneModel` subclass via :meth:`GeneModel.from_gene`, which requires the
gene's CDS coordinates plus a reference sequence and, from them, derives the
coding sequence and four sequence attributes (:attr:`~GeneModel.protein`,
:attr:`~GeneModel.rna`, :attr:`~GeneModel.dna` and
:attr:`~GeneModel.revcomp_dna`)."""

from theiagene.lib.sequence import complement, reverse_complement, translate


# attribute keys under which CDS coordinate segments are filed (case variants)
CDS_KEYS = ("CDS", "cds")


class Gene:
    """Coordinates and metadata for a single query gene.

    Coordinates are 0-based, half-open and refer to the reference contig.
    ``parts`` is a ``{feature_type: [(start, end), ...]}`` dict -- CDS segments
    live under a ``"CDS"`` key (see :attr:`cds`); the derived properties order
    those into translation (5'->3') order using ``strand`` (1/-1/None).
    ``qualifiers`` holds the parsed feature attributes ({name: [values]}) used to
    match the gene against a query set, ``contig_candidates`` the seqid/name(s) a
    resolved ``contig`` is chosen from, and ``contig_seq`` the reference sequence
    when one is available (required to build a :class:`GeneModel`)."""

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
    ):
        self.gene_id = gene_id
        self.contig = contig
        self.strand = strand
        self.transl_table = transl_table
        self.product = product if product is not None else gene_id
        # parts must be a {feature_type: [(start, end), ...]} dict (0-based,
        # half-open); CDS segments live under a "CDS" key (see the cds property)
        self.parts = {}
        if parts:
            if not isinstance(parts, dict):
                raise TypeError(
                    "Gene 'parts' must be a {feature_type: [(start, end), ...]} "
                    f"dict, not {type(parts).__name__}"
                )
            for feature, segments in parts.items():
                self.parts[feature] = [(int(s), int(e)) for s, e in segments]
        self.qualifiers = qualifiers if qualifiers is not None else {}
        self.contig_seq = contig_seq
        self.contig_candidates = (
            list(contig_candidates)
            if contig_candidates is not None
            else ([contig] if contig is not None else [])
        )

    @property
    def cds(self) -> list:
        """CDS coordinate segments (case-insensitive key lookup); ``[]`` if none"""
        for key in CDS_KEYS:
            if key in self.parts:
                return self.parts[key]
        return []

    def add_part(self, start, end, feature: str = "CDS") -> None:
        """Append a genomic ``(start, end)`` segment (0-based, half-open) under
        ``feature`` (the CDS coordinate list by default)"""
        self.parts.setdefault(feature, []).append((int(start), int(end)))

    @property
    def genomic_start(self):
        """Left-most genomic coordinate across the CDS parts (None if empty)"""
        return min((s for s, _ in self.cds), default=None)

    @property
    def genomic_end(self):
        """Right-most genomic coordinate across the CDS parts (None if empty)"""
        return max((e for _, e in self.cds), default=None)

    @property
    def genomic_positions(self):
        """Every coding-base genomic position in translation (5'->3') order"""
        positions = []
        for start, end in sorted(self.cds):
            positions.extend(range(start, end))
        if self.strand == -1:
            positions.reverse()
        return positions


class GeneModel(Gene):
    """A :class:`Gene` plus its reference sequence and the derived sequences.

    :meth:`finalize` translates the gene's coordinates against a contig sequence
    into the coding sequence (``ref_coding``, in 5'->3' translation order,
    reverse-complemented for minus-strand genes) and its ``pos2cds`` index --
    both used by the variant placement/codon math in ``theiagene.lib.variant`` --
    together with four sequence attributes:

    ``protein``      translation of the spliced coding sequence
    ``rna``          the spliced coding sequence as RNA (coding strand, T->U)
    ``dna``          the full gene span, introns included, on the coding strand
    ``revcomp_dna``  the reverse complement of ``dna`` (template strand)
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ref_coding = ""
        self.pos2cds = {}
        self.protein = ""
        self.rna = ""
        self.dna = ""
        self.revcomp_dna = ""

    @classmethod
    def from_gene(cls, gene: Gene) -> "GeneModel":
        """Generate a finalized :class:`GeneModel` from a :class:`Gene`.

        The model is derived entirely from the gene: it requires the gene to
        carry CDS coordinates (a ``CDS``/``cds`` :attr:`~Gene.cds` list) and a
        reference :attr:`~Gene.contig_seq`, raising ``ValueError`` if either is
        missing.  Callers resolve the gene's identity (``gene_id``, ``product``,
        ``contig``, ``transl_table``) on the gene itself and let the model copy
        it, so there is a single place that decides those values."""
        if not gene.cds:
            raise ValueError(f"cannot model '{gene.gene_id}': no CDS coordinates")
        if gene.contig_seq is None:
            raise ValueError(f"cannot model '{gene.gene_id}': no reference sequence")
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

    def finalize(self, contig_seq: str = None) -> None:
        """Build the coding sequence, position index and derived sequences.

        Uses ``contig_seq`` when given, otherwise the gene's stored
        :attr:`~Gene.contig_seq`; raises ``ValueError`` when neither is set."""
        if contig_seq is None:
            contig_seq = self.contig_seq
        if contig_seq is None:
            raise ValueError(f"cannot finalize '{self.gene_id}': no reference sequence")
        self.contig_seq = contig_seq
        positions = self.genomic_positions
        coding = "".join(contig_seq[p].upper() for p in positions)
        if self.strand == -1:
            # positions are already reversed; complement each base to get revcomp
            coding = complement(coding)
        self.ref_coding = coding
        self.pos2cds = {g: i for i, g in enumerate(positions)}
        self.protein = translate(coding, self.transl_table)
        self.rna = coding.replace("T", "U")
        # full gene span including introns, on the coding strand
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
