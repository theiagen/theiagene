# theiagene

theiagene is a gene-centric data manipulation toolkit and library. 
Provides a single `theiagene` command-line entrypoint with two
subcommands:

- **`gene_coverage`** — quantify the breadth and depth of coverage over query
  genes from a BAM
- **`variant_annotation`** — annotate the protein-level consequences
  (missense / synonymous / nonsense substitutions, in-frame insertions/deletions
  and frameshifts) of variants that overlap query genes

## Installation

```bash
pip install .
# or, for development:
pip install -e '.[test]'
```

## Usage

```bash
theiagene --help
theiagene gene_coverage --help
theiagene variant_annotation --help
```

### gene_coverage

Report average depth and percent coverage per query gene. Coordinates come from
a reference GenBank/GFF or a BED file; outputs are written to the working
directory as `DEPTH_DICT.json`, `COVERAGE_DICT.json` and `COVERAGE_STATS.tsv`.

```bash
theiagene gene_coverage \
  --bam sample.sorted.bam \
  --bedfile regions.bed

theiagene gene_coverage \
  --bam sample.sorted.bam \
  --reference_gbff reference.gbff \
  --query_genes FKS1 ERG11
```

### variant_annotation

Annotate the protein-level effect of each gene-overlapping variant. Accepts a
raw VCF, together with a reference `--reference_gbff` (or `--reference_gff` and `--reference_fa`). Extracts gene coordinates from VCF into a `GENE_VARIANTS.vcf` and writes
the report to `VARIANT_ANNOTATIONS.txt` (override with `--output`).

```bash
theiagene variant_annotation \
  --vcf GENE_VARIANTS.vcf \
  --reference_gbff reference.gbff \
  --query_genes "lanosterol 14-alpha demethylase" \
  --exact_match
```

Example report line:

```
lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.395A>T p.Tyr132Phe; A:20 T:0)
```

NOTE: frameshift mutations will terminate variant calls that are 3' relative to the frameshift.

## HGVS notation

Every annotation carries a nucleotide description (`c.`) and a protein
description (`p.`) in [HGVS](https://hgvs-nomenclature.org/) form. Both are
1-based and relative to the gene rather than the contig: `c.1` is the first base
of the start codon and `p.1` the initiator methionine, and for minus-strand
genes bases are reported on the coding strand (complemented, 5'→3'). Amino acids
use three-letter codes, with `Ter` for a stop codon.

`c.` — position in the coding sequence, then the change:

| Example | Meaning |
| --- | --- |
| `c.395A>T` | A at coding base 395 substituted by T |
| `c.395_396insTTC` | TTC inserted between bases 395 and 396 |
| `c.395delA` | the single base A at 395 deleted |
| `c.395_397del` | bases 395–397 deleted |
| `c.395delinsTT` | base 395 replaced by TT |
| `c.395_397delinsTT` | bases 395–397 replaced by TT |

`p.` — affected residue(s), then the consequence:

| Example | Meaning | SO term |
| --- | --- | --- |
| `p.Tyr132Phe` | Tyr132 becomes Phe | `missense_variant` |
| `p.Tyr132=` / `p.=` | protein product unchanged | `synonymous_variant` |
| `p.Tyr132Ter` | Tyr132 becomes a stop codon | `stop_gained` |
| `p.Ter132Tyr` / `p.Ter132ext` | stop lost, translation reads through | `stop_lost` |
| `p.Met1Thr` | start codon lost | `start_lost` |
| `p.Tyr132fs` | reading frame shifts from Tyr132 onward | `frameshift_variant` |
| `p.Tyr132del` / `p.Tyr132_Phe134del` | one / several residues removed | `inframe_deletion` |
| `p.Tyr132_Phe133insLeuVal` | LeuVal inserted between Tyr132 and Phe133 | `inframe_insertion` |
| `p.(=)` | change lies past the stop codon; product unchanged | `stop_retained_variant` |
| `p.?` | consequence undetermined (truncated final codon) | `coding_sequence_variant` |

Only coding changes are described, so there is no `g.`/`n.` notation and no
intronic offsets (`c.123+4`); a variant touching no coding base is skipped.
Indels are left-aligned as they appear in the VCF rather than shifted to the
3'-most position HGVS prescribes, and a duplication is reported as an insertion
rather than `dup`.