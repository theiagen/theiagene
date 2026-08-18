# theiagene

theiagene is a gene-centric data manipulation toolkit and library. 
Provides a single `theiagene` command-line entrypoint with three
subcommands:

- **`gene_coverage`** — quantify the breadth, depth, and read support of
  coverage over query genes from a BAM
- **`extract_variants`** — extract a sub-VCF of variants that fall within query
  genes from a VCF
- **`report_variants`** — render VEP variant annotations into product-named
  report lines

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
theiagene extract_variants --help
theiagene report_variants --help
```

### gene_coverage

Report average depth, percent coverage, and mapped reads per query gene.
Coordinates come from a reference GFF or a BED file; outputs are written to the
working directory as `DEPTH_DICT.json`, `COVERAGE_DICT.json`, `READS_DICT.json`,
`READS_PASS_DICT.json` and `COVERAGE_STATS.tsv`. A gene whose mapped reads fall
below `--min_reads_mapped` (default 1) is flagged as failing in
`READS_PASS_DICT.json`, but keeps its measured depth, breadth, and read count —
the flag is a call layered onto the measurements, not a filter of them.

All three measurements are taken over one read set, so `reads_mapped` describes
the alignments the depth and breadth were computed from:

| filter | effect |
| --- | --- |
| always applied | unmapped, secondary, supplementary, QC_fail and duplicate alignments are excluded |
| `--min_mapping_quality` (default 0) | an alignment below it counts toward neither reads nor depth/breadth |
| `--min_base_quality` (default 0) | a base below it counts toward neither depth nor breadth, and a read with no qualifying base in a region is not mapped to it |
| `--min_depth` (default 1) | a base at or above this depth counts as covered, setting breadth |

`--min_base_quality` was previously `--min_quality`, which still works as a
deprecated alias.

```bash
theiagene gene_coverage \
  --bam sample.sorted.bam \
  --bedfile regions.bed

theiagene gene_coverage \
  --bam sample.sorted.bam \
  --reference_gff reference.gff \
  --query_genes FKS1 ERG11
```

### extract_variants

Write a sub-VCF containing only the variants that overlap the `feature_type`
(CDS by default) segments of the query genes. Coordinates come from a reference
GFF or a BED file; each kept record is annotated with the overlapping query
name(s) in a `GENE` INFO field. Output defaults to `EXTRACTED_VARIANTS.vcf`.

```bash
theiagene extract_variants \
  --vcf sample.vcf \
  --reference_gff reference.gff \
  --query_genes FKS1 ERG11

theiagene extract_variants \
  --vcf sample.vcf \
  --bedfile regions.bed
```

### report_variants

Render a VEP `--tab` output TSV into gene-labelled report lines. Rows with a
suppressed consequence, no HGVSc/HGVSp string, or a feature that resolves to no
CDS product are dropped. Each remaining row becomes a gene label, the quoted CDS
product resolved through the reference GFF, and the consequence with the
transcript/protein prefixes stripped from its HGVS strings, e.g.:

```
ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg)
```

The label is the `--query_genes` term (or, failing that, the `--bedfile` name)
matching the row's feature — the name that was asked about rather than the full
product it resolved to. With no query list, or when none of its terms match, the
label falls back to the normalized product name
(`lanosterol.14-alpha.demethylase: "lanosterol 14-alpha demethylase" (...)`).
The product is quoted so a name carrying commas survives being joined into a
comma-delimited report.

Lines print to stdout unless `--output` is given. When the source VCF is passed
via `--vcf`, each line also carries the variant's per-allele read depths (e.g.
`; T:0 C:562`).

```bash
theiagene report_variants \
  --vep_tsv variants.vep.tsv \
  --reference_gff reference.gff

theiagene report_variants \
  --vep_tsv variants.vep.tsv \
  --reference_gff reference.gff \
  --vcf sample.vcf \
  --suppress synonymous_variant \
  --output VARIANT_REPORT.txt
```

## Library

The subcommands share a gene/feature data model — the `Feature` and
`FeatureCol` classes — that turns a flat GFF3 annotation into a
navigable gene → RNA → CDS/exon hierarchy. See
[src/theiagene/lib/README.md](src/theiagene/lib/README.md) for a human-readable
introduction and full API reference.