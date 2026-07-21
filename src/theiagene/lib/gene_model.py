"""Gene coordinate model (:class:`Gene`) and its sequence-bearing subclass.

A :class:`Gene` is the shared coordinate currency of theiagene: it captures the
metadata and genomic segments of one query gene, without any reference sequence.
The coordinate parsers in ``theiagene.lib.parsers`` produce one per matched
feature -- ``gene_coverage`` uses these directly (it only needs coordinates),
while ``variant_annotation`` produces the :class:`GeneModel` subclass, which
attaches the reference sequence and, from it, derives the coding sequence and
four sequence attributes (:attr:`~GeneModel.protein`, :attr:`~GeneModel.rna`,
:attr:`~GeneModel.dna` and :attr:`~GeneModel.revcomp_dna`)."""

from theiagene.lib.sequence import complement, reverse_complement, translate


class Gene:
    """Coordinates and metadata for a single query gene, without sequence.

    Coordinates are 0-based, half-open and refer to the reference contig.
    ``parts`` holds one ``(start, end)`` segment per exon/CDS piece; the derived
    properties order them into translation (5'->3') order using ``strand``
    (1/-1/None)."""

    def __init__(
        self,
        gene_id,
        contig,
        strand=None,
        transl_table=1,
        product=None,
        parts=None,
    ):
        self.gene_id = gene_id
        self.contig = contig
        self.strand = strand
        self.transl_table = transl_table
        self.product = product if product is not None else gene_id
        self.parts = [(int(s), int(e)) for s, e in (parts or [])]

    def add_part(self, start, end) -> None:
        """Append a genomic ``(start, end)`` segment (0-based, half-open)"""
        self.parts.append((int(start), int(end)))

    @property
    def genomic_start(self):
        """Left-most genomic coordinate across all parts (None if empty)"""
        return min((s for s, _ in self.parts), default=None)

    @property
    def genomic_end(self):
        """Right-most genomic coordinate across all parts (None if empty)"""
        return max((e for _, e in self.parts), default=None)

    @property
    def genomic_positions(self):
        """Every coding-base genomic position in translation (5'->3') order"""
        positions = []
        for start, end in sorted(self.parts):
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

    def finalize(self, contig_seq: str) -> None:
        """Build the coding sequence, position index and derived sequences"""
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
