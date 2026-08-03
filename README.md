# theiagene

theiagene is a gene-centric data manipulation toolkit and library. 
Provides a single `theiagene` command-line entrypoint with three
subcommands:

- **`gene_coverage`** — quantify the breadth and depth of coverage over query
  genes from a BAM
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

Render a VEP `--tab` output TSV into product-named report lines. Rows with a
suppressed consequence, no HGVSc/HGVSp string, or an unresolvable feature are
dropped; the remaining rows have their transcript/protein prefixes rewritten to
the CDS product name resolved through the reference GFF, e.g.:

```
lanosterol.14-alpha.demethylase: lanosterol 14-alpha demethylase (missense_variant c.428A>G p.Lys143Arg)
```

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
`FeatureCol` classes — that turns a flat GFF/GenBank annotation into a
navigable gene → RNA → CDS/exon hierarchy. See
[src/theiagene/lib/README.md](src/theiagene/lib/README.md) for a human-readable
introduction and full API reference.