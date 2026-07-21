"""Per-gene coding-sequence model shared by the variant-annotation backend.

A :class:`GeneModel` captures the coding sequence, strand, genomic span and
translation table of one query gene, together with the position index and
reference protein needed to place and translate variants against it.  The GBFF
and GFF model builders in ``theiagene.variant_annotation`` populate one per
query gene; :class:`theiagene.lib.variant.Variant` annotates changes against
it."""

from theiagene.lib.sequence import complement, translate


class GeneModel:
    """A coding sequence model for a single query gene.

    Coordinates are 0-based, half-open and refer to the reference contig.
    ``genomic_positions`` lists every coding base position in translation
    (5'->3') order; ``ref_coding`` is the reference coding sequence in that same
    order (reverse-complemented for minus-strand genes); ``pos2cds`` maps a
    genomic position to its index within ``ref_coding``."""

    def __init__(self, gene_id, product, contig, strand, transl_table):
        self.gene_id = gene_id
        self.product = product
        self.contig = contig
        self.strand = strand
        self.transl_table = transl_table
        self.genomic_positions = []
        self.pos2cds = {}
        self.ref_coding = ""
        self.ref_protein = ""
        # genomic span used for interval-overlap gene assignment
        self.genomic_start = None
        self.genomic_end = None

    def finalize(self, contig_seq: str) -> None:
        """Build the coding sequence, position index and reference protein"""
        bases = [contig_seq[p].upper() for p in self.genomic_positions]
        coding = "".join(bases)
        if self.strand == -1:
            # positions are already reversed; complement each base to get revcomp
            coding = complement(coding)
        self.ref_coding = coding
        self.pos2cds = {g: i for i, g in enumerate(self.genomic_positions)}
        self.ref_protein = translate(coding, self.transl_table)
        if self.genomic_positions:
            self.genomic_start = min(self.genomic_positions)
            self.genomic_end = max(self.genomic_positions) + 1

    def codon(self, codon_number: int) -> str:
        """Return the reference codon (1-based codon number) or '' if incomplete"""
        start = (codon_number - 1) * 3
        codon = self.ref_coding[start : start + 3]
        return codon if len(codon) == 3 else ""

    def aa_at(self, codon_number: int) -> str:
        """Return the reference amino acid (one letter) at a 1-based codon"""
        codon = self.codon(codon_number)
        return translate(codon, self.transl_table) if codon else "X"
