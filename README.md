# theiagene

theiagene is a gene-centric data manipulation toolkit and library. 
Provides a single `theiagene` command-line entrypoint with two
subcommands:

- **`gene_coverage`** — quantify the breadth and depth of coverage over query
  genes from a BAM

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