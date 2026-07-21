"""Nucleotide/protein sequence helpers shared by the theiagene commands."""

from Bio.Seq import Seq
from Bio.Data import IUPACData


# single-letter -> three-letter amino acid codes, with HGVS extensions
_AA_1TO3 = {k.upper(): v for k, v in IUPACData.protein_letters_1to3.items()}
_AA_1TO3["*"] = "Ter"
_AA_1TO3["X"] = "Xaa"

_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

_VALID_BASES = set("ACGTNacgtn")


def is_nucleotide_allele(allele) -> bool:
    """True only for plain nucleotide alleles.

    Filters symbolic alleles ('<DEL>', '<*>'), the spanning-deletion allele
    ('*') emitted by GATK/bcftools, and empty/missing alleles, all of which
    must not be translated as substitutions."""
    return bool(allele) and all(base in _VALID_BASES for base in allele)


def aa3(aa: str) -> str:
    """Return the three-letter code for a one-letter amino acid symbol"""
    return _AA_1TO3.get(aa.upper(), "Xaa")


def complement(seq: str) -> str:
    """Complement (not reverse) a nucleotide string"""
    return seq.translate(_COMPLEMENT)


def translate(seq: str, table) -> str:
    """Translate a nucleotide string, truncating to a whole number of codons"""
    trimmed = seq[: len(seq) - (len(seq) % 3)]
    if not trimmed:
        return ""
    return str(Seq(trimmed).translate(table=table))
