"""Protein-level annotation of a variant against a gene's coding model.

A :class:`Variant` pairs a normalized REF/ALT change with the
:class:`theiagene.lib.gene_model.GeneModel` it overlaps and classifies the
protein consequence (Sequence Ontology term plus HGVS c./p. notation).  The
module-level helpers ``_protein_consequence``, ``_first_diff`` and ``_res`` are
the pure protein-comparison routines it relies on."""

import logging

from theiagene.lib.sequence import aa3, complement, translate
from theiagene.lib.gene_model import GeneModel


logger = logging.getLogger(__name__)


class Variant:
    """A single normalized nucleotide change evaluated against one gene.

    A ``Variant`` pairs a trimmed REF/ALT change -- ``ref_seg``/``alt_seg`` at
    0-based genomic ``changed_pos0``, as produced by ``normalize_indel`` -- with
    the :class:`GeneModel` it is annotated against (``model``).  Its methods are
    the transformations that classify the change and build its HGVS c./p.
    notation.

    ``annotate`` returns an annotation dict with keys ``so_term``, ``hgvs_c``,
    ``hgvs_p``, ``cds_pos`` and ``is_frameshift``, or ``None`` when the change
    does not fall on the model's coding sequence."""

    def __init__(
        self, model: GeneModel, changed_pos0: int, ref_seg: str, alt_seg: str
    ):
        self.model = model
        # position of variant
        self.changed_pos0 = changed_pos0
        self.ref_seg = ref_seg
        self.alt_seg = alt_seg

    def annotate(self) -> dict:
        """Classify the change against the model.

        Dispatches to :meth:`annotate_snp` for a single-base substitution and to
        :meth:`annotate_indel` otherwise; returns ``None`` for a no-op change
        (empty REF and ALT segments) or one that does not touch coding bases."""
        if len(self.ref_seg) == 1 and len(self.alt_seg) == 1:
            return self.annotate_snp()
        # if no SNPs are annotated by no ref and alt
        if not self.ref_seg and not self.alt_seg:
            return None
        # guard against indel path for erroneously labeled SNVs that are not variants 
        elif len(self.ref_seg) > 1 and self.ref_seg == self.alt_seg:
            return None
        return self.annotate_indel()

    def annotate_snp(self) -> dict:
        """Annotate a single-nucleotide substitution against the coding sequence"""
        model = self.model
        pos0, ref, alt = self.changed_pos0, self.ref_seg, self.alt_seg
        cds_idx = model.pos2cds.get(pos0)
        if cds_idx is None:
            return None  # not within a coding base (e.g. intronic anchor)

        # account for reverse complement SNPs/SNVs
        coding_ref = ref if model.strand == 1 else complement(ref)
        coding_alt = alt if model.strand == 1 else complement(alt)

        expected = model.ref_coding[cds_idx]
        if expected != coding_ref:
            logger.warning(
                f"{model.gene_id}: reference base '{expected}' at c.{cds_idx + 1} "
                f"disagrees with VCF REF '{coding_ref}' (coding strand)"
            )

        codon_number = cds_idx // 3 + 1
        pos_in_codon = cds_idx % 3
        ref_codon = model.codon(codon_number)
        hgvs_c = f"c.{cds_idx + 1}{coding_ref}>{coding_alt}"

        if not ref_codon:
            # incomplete terminal codon; report nucleotide change only
            return {
                "so_term": "coding_sequence_variant",
                "hgvs_c": hgvs_c,
                "hgvs_p": "p.?",
                "cds_pos": cds_idx,
                "is_frameshift": False,
            }

        mut_codon = ( 
            ref_codon[:pos_in_codon] + coding_alt + ref_codon[pos_in_codon + 1 :]
        )
        ref_aa = translate(ref_codon, model.transl_table)
        alt_aa = translate(mut_codon, model.transl_table)

        if ref_aa == alt_aa:
            so_term = "synonymous_variant"
            hgvs_p = f"p.{aa3(ref_aa)}{codon_number}="
        elif alt_aa == "*":
            so_term = "stop_gained"
            hgvs_p = f"p.{aa3(ref_aa)}{codon_number}Ter"
        elif ref_aa == "*":
            so_term = "stop_lost"
            hgvs_p = f"p.Ter{codon_number}{aa3(alt_aa)}"
        elif codon_number == 1 and ref_aa == "M":
            so_term = "start_lost"
            hgvs_p = f"p.{aa3(ref_aa)}1{aa3(alt_aa)}"
        else:
            so_term = "missense_variant"
            hgvs_p = f"p.{aa3(ref_aa)}{codon_number}{aa3(alt_aa)}"

        return {
            "so_term": so_term,
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "cds_pos": cds_idx,
            "is_frameshift": False,
        }

    def annotate_indel(self) -> dict:
        """Annotate an insertion, deletion or delins against the coding sequence"""
        model = self.model
        changed_pos0, ref_seg, alt_seg = (
            self.changed_pos0,
            self.ref_seg,
            self.alt_seg,
        )
        # coding indices of the deleted reference bases (those that are coding)
        deleted_idx = sorted(
            model.pos2cds[changed_pos0 + i]
            for i in range(len(ref_seg))
            if (changed_pos0 + i) in model.pos2cds
        )
        coding_alt = alt_seg if model.strand == 1 else complement(alt_seg)
        if model.strand == -1:
            coding_alt = coding_alt[::-1]  # reverse to restore 5'->3' coding order

        deleted_cds = len(deleted_idx)
        inserted_cds = len(coding_alt)
        frame_delta = (inserted_cds - deleted_cds) % 3

        # locate the variant within the coding sequence for splicing / reporting
        if deleted_idx:
            cds_start = deleted_idx[0]
            cds_end = deleted_idx[-1] + 1
        else:
            # pure insertion: require a coding base immediately on one side, then
            # splice at the count of coding bases lying 5' of the insertion point.
            # This is strand-aware and correct at exon/CDS boundaries
            if (
                changed_pos0 - 1
            ) not in model.pos2cds and changed_pos0 not in model.pos2cds:
                return None  # insertion is not adjacent to any coding base
            if model.strand == 1:
                cds_start = sum(1 for g in model.cds_positions if g < changed_pos0)
            else:
                cds_start = sum(1 for g in model.cds_positions if g >= changed_pos0)
            cds_end = cds_start

        first_codon = cds_start // 3 + 1

        if frame_delta != 0:
            return {
                "so_term": "frameshift_variant",
                "hgvs_c": self._hgvs_c_indel(cds_start, cds_end, coding_alt),
                "hgvs_p": f"p.{aa3(model.aa_at(first_codon))}{first_codon}fs",
                "cds_pos": cds_start,
                "is_frameshift": True,
            }

        # in-frame: build the mutant coding sequence and compare proteins
        mut_coding = (
            model.ref_coding[:cds_start] + coding_alt + model.ref_coding[cds_end:]
        )
        mut_protein = translate(mut_coding, model.transl_table)
        ref_protein = model.protein

        hgvs_c = self._hgvs_c_indel(cds_start, cds_end, coding_alt)
        so_term, hgvs_p = _protein_consequence(
            ref_protein, mut_protein, first_codon, deleted_cds, inserted_cds
        )
        return {
            "so_term": so_term,
            "hgvs_c": hgvs_c,
            "hgvs_p": hgvs_p,
            "cds_pos": cds_start,
            "is_frameshift": False,
        }

    def _hgvs_c_indel(self, cds_start, cds_end, coding_alt) -> str:
        """Build the HGVS c. string for an indel/delins in coding coordinates"""
        model = self.model
        ref_seg = self.ref_seg
        deleted_len = cds_end - cds_start
        if not ref_seg or deleted_len == 0:
            # pure insertion between cds_start-1 and cds_start (1-based cds_start, +1)
            return f"c.{cds_start}_{cds_start + 1}ins{coding_alt}"
        deleted_bases = model.ref_coding[cds_start:cds_end]
        if not coding_alt:
            if deleted_len == 1:
                return f"c.{cds_start + 1}del{deleted_bases}"
            return f"c.{cds_start + 1}_{cds_end}del"
        # delins
        if deleted_len == 1:
            return f"c.{cds_start + 1}delins{coding_alt}"
        return f"c.{cds_start + 1}_{cds_end}delins{coding_alt}"


def _protein_consequence(
    ref_protein, mut_protein, first_codon, deleted_cds, inserted_cds
) -> tuple:
    """Classify an in-frame protein change and build the HGVS p. string"""
    net = inserted_cds - deleted_cds
    net_aa = net // 3  # frame_delta == 0, so net is a whole number of codons

    # translated products up to (and excluding) the first stop codon
    ref_prod = ref_protein.split("*", 1)[0]
    mut_prod = mut_protein.split("*", 1)[0]
    first_res = first_codon - 1  # 0-based first affected residue

    # a variant lying entirely 3' of the reference stop codon cannot alter the
    # protein product; only report a change if it actually reaches the terminator
    if first_res >= len(ref_prod):
        if mut_prod == ref_prod:
            return "stop_retained_variant", "p.(=)"
        return "stop_lost", f"p.Ter{len(ref_prod) + 1}ext"

    # a genuine stop_gained truncates the product earlier than the indel's own
    # length change accounts for; a plain in-frame deletion only shortens it by
    # net_aa residues and must not be mistaken for a new stop. Anchor the Ter at
    # the mutant's first stop, not at the first *changed* residue (which may be a
    # sense codon upstream of the new stop in a multi-codon delins)
    if len(mut_prod) < len(ref_prod) + net_aa:
        stop_pos = mut_protein.find("*")
        codon = stop_pos + 1
        ref_res = _res(ref_protein, stop_pos - net_aa)
        return "stop_gained", f"p.{aa3(ref_res)}{codon}Ter"

    idx = _first_diff(ref_protein, mut_protein)

    if net < 0:  # net deletion of amino acids
        n_del = (-net) // 3
        first = idx
        last = idx + n_del - 1
        if n_del <= 1:
            hgvs_p = f"p.{aa3(_res(ref_protein, first))}{first + 1}del"
        else:
            hgvs_p = (
                f"p.{aa3(_res(ref_protein, first))}{first + 1}_"
                f"{aa3(_res(ref_protein, last))}{last + 1}del"
            )
        return "inframe_deletion", hgvs_p

    if net > 0:  # net insertion of amino acids
        n_ins = net // 3
        before = idx - 1
        after = idx
        inserted = "".join(aa3(_res(mut_protein, idx + k)) for k in range(n_ins))
        hgvs_p = (
            f"p.{aa3(_res(ref_protein, before))}{before + 1}_"
            f"{aa3(_res(ref_protein, after))}{after + 1}ins{inserted}"
        )
        return "inframe_insertion", hgvs_p

    # equal length: one or more codons substituted (MNV / in-frame delins)
    if idx >= len(ref_protein) or idx >= len(mut_protein):
        return "synonymous_variant", "p.="
    if _res(ref_protein, idx) == _res(mut_protein, idx):
        return "synonymous_variant", "p.="
    return (
        "missense_variant",
        f"p.{aa3(_res(ref_protein, idx))}{idx + 1}{aa3(_res(mut_protein, idx))}",
    )


def _first_diff(a: str, b: str) -> int:
    """Index of the first position at which two strings differ"""
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return i
    return min(len(a), len(b))


def _res(protein: str, idx: int) -> str:
    """Residue at ``idx`` (one letter), or 'X' if out of range"""
    return protein[idx] if 0 <= idx < len(protein) else "X"
