# theiagene

Gene coverage and variant annotation toolkit for Theiagen bioinformatics
workflows. Provides a single `theiagene` command-line entrypoint with two
subcommands:

- **`gene_coverage`** — quantify the breadth and depth of coverage over query
  genes from a BAM (a dependency of
  [`task_gene_coverage.wdl`](https://github.com/theiagen/public_health_bioinformatics)).
- **`variant_annotation`** — annotate the protein-level consequences
  (missense / synonymous / nonsense substitutions, in-frame insertions/deletions
  and frameshifts) of variants that overlap query genes.

These were previously two standalone scripts (`gene_coverage.py`,
`variant_annotation.py`) shipped inside the `pysam` Docker build; they are now a
proper package with the functionality they shared factored into libraries.

## Installation

```bash
pip install .
# or, for development:
pip install -e '.[test]'
```

Runtime dependencies are `pysam` and `biopython` (Python ≥ 3.9).

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
Supplying `--vcf` additionally extracts gene-overlapping variants to
`GENE_VARIANTS.vcf`.

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
raw VCF or a `GENE_VARIANTS.vcf` produced by `gene_coverage`, together with a
reference `--reference_gbff` (or `--reference_gff` and `--reference_fa`). Writes
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

The two commands can be run without installing via
`python -m theiagene <command> ...`.

## Layout

```
src/theiagene/
├── cli.py                 # unified "theiagene" entrypoint (subcommand dispatch)
├── gene_coverage.py       # gene_coverage subcommand
├── variant_annotation.py  # variant_annotation subcommand
└── lib/                   # functionality shared by the two commands
    ├── logging_config.py  # shared logging setup
    ├── query.py           # query-gene identifier matching / normalization
    ├── sequence.py        # nucleotide/protein sequence helpers
    ├── vcf.py             # VCF coordinate flattening + gene extraction
    └── io_utils.py        # WDL-compatible JSON output
tests/                     # pytest suite (mirrors src, plus CLI + lib tests)
```

Each subcommand module exposes `add_arguments(parser)` and `run_cli(args)`, which
`cli.py` wires together; the shared helpers each command uses are imported from
`theiagene.lib` (and re-exported on the command modules so existing call sites
keep working).

## Testing

```bash
pytest
```

The suite has no network dependency; VCF/GBFF/BAM fixtures are synthesized in
`tmp_path`.
