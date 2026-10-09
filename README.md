# theiagene

theiagene is a gene-centric data manipulation toolkit and library. 
Provides a single `theiagene` command-line entrypoint with four
subcommands:

- **`gene_coverage`** — quantify the breadth, depth, and read support of
  coverage over query genes from a BAM
- **`extract_variants`** — extract a sub-VCF of variants that fall within query
  genes from a VCF
- **`report_variants`** — render SnpEff variant annotations into a product-named
  report TSV
- **`prepare_snpeff`** — prepare a reference GFF and config for building a
  SnpEff database

## Installation

```bash
pip install .
# or, for development:
pip install -e '.[test]'
```

## Usage

### Query and coordinate sources

Query and coordinate arguments are handled hierarchically:

| given | query source | coordinate source |
| --- | --- | --- |
| `--query_genes` + `--reference_gff` | query terms | GFF |
| `--query_genes` + `--reference_gff` + `--bedfile` | query terms | GFF — BED is ignored |
| `--query_genes` + `--bedfile` | query terms | corresponding BED coordinate columns |
| `--reference_gff` + `--bedfile` | BED name column | GFF |
| `--bedfile` alone | BED name column | BED coordinate columns |

A GFF is the preferred query coordinate source, and `--query_genes` is the preferred
query name source — a BED steps in for whichever may be missing. At least one
coordinate source and one query source is required (`--bedfile` alone satisfies
both); `report_variants` _requires_ a GFF, so its `--bedfile` only
ever supplies names.

GFF query matching is case-insensitive and accommodates substrings (unless
`--exact_match` is specified), whereas BED rows are selected by an exact, 
case-sensitive match on the name column. `--query_genes erg11` therefore 
finds `ERG11` in a GFF but not in a BED.


### Subcommands

```bash
theiagene --help
theiagene gene_coverage --help
theiagene extract_variants --help
theiagene report_variants --help
theiagene prepare_snpeff --help
```

#### gene_coverage

Report average depth, percent coverage, mapped reads, and quantified length per
query gene, over the coordinates resolved as described above; outputs are
written to the working directory as `DEPTH_DICT.json`, `COVERAGE_DICT.json`,
`READS_DICT.json`, `READS_PASS_DICT.json`, `LENGTHS_DICT.json` and
`COVERAGE_STATS.tsv`. A gene whose mapped reads fall below `--min_reads_mapped`
(default 1) is flagged as failing in `READS_PASS_DICT.json`, but keeps its
measured depth, breadth, and read count.

`LENGTHS_DICT.json` gives the number of reference bases each gene's depth and breadth
were computed over — the denominator of both ratios. Its bases are counted after
the gene's ranges are merged, so a base shared by two overlapping segments (e.g.
isoform CDS) contributes once, and a gene split across several segments or
contigs reports their combined length.

A query gene that resolves to no coordinates in the annotation is reported as
`NA` in all five outputs.

Each measurement is also summarized across genes into a single-value file:
`MEAN_DEPTH`, `MEAN_COVERAGE`, `MEAN_READS`, `TOTAL_DEPTH`, `TOTAL_COVERAGE` and
`TOTAL_READS`. Depth and breadth are per-base quantities, so their means weight
each gene by its quantified length — the mean over every quantified base rather
than the mean of per-gene values; reads are counted per gene, so their mean
counts each gene once. A gene reported as `NA` was never measured, so it enters
neither the mean nor the total; when nothing was measured at all, both files are
blank rather than `0`.

| filter | effect |
| --- | --- |
| always applied | unmapped, secondary, supplementary, QC_fail and duplicate alignments are excluded |
| `--min_mapping_quality` (default 0) | an alignment below it counts toward neither reads nor depth/breadth |
| `--min_base_quality` (default 0) | a base below it counts toward neither depth nor breadth, and a read with no qualifying base in a region is not mapped to it |
| `--min_depth` (default 1) | a base at or above this depth counts as covered, setting breadth |

```bash
theiagene gene_coverage \
  --bam sample.sorted.bam \
  --bedfile regions.bed

theiagene gene_coverage \
  --bam sample.sorted.bam \
  --reference_gff reference.gff \
  --query_genes FKS1,ERG11
```

#### extract_variants

Write a sub-VCF containing only the variants that overlap the `feature_type`
(CDS by default) segments of the query genes, over the coordinates resolved as
described above. Each kept record is annotated with the query that retrieved it
in a `GENE` INFO field. Output defaults to `EXTRACTED_VARIANTS.vcf`.

```bash
theiagene extract_variants \
  --vcf sample.vcf \
  --reference_gff reference.gff \
  --query_genes FKS1,ERG11

theiagene extract_variants \
  --vcf sample.vcf \
  --bedfile regions.bed
```

#### report_variants

Render a SnpEff-annotated VCF into a gene-labelled report TSV, one record per
kept `ANN` entry. `--reference_gff` should be the GFF the SnpEff database was
built from (see [prepare_snpeff](#prepare_snpeff)), so the transcripts SnpEff
names are the ones resolved here:

| column | content | example |
| --- | --- | --- |
| `GENE` | gene label (see below) | `ERG11` |
| `HGVSc` | SnpEff's HGVS.c, prefixed with the transcript ID | `rna-x:c.428A>G` |
| `HGVSp` | SnpEff's HGVS.p, prefixed with the CDS `protein_id` (else the transcript ID) | `prot-x:p.Lys143Arg` |
| `NT` | abbreviated nucleotide change: HGVS form without its prefix | `428A>G`, `383delA` |
| `AA` | abbreviated amino acid change: one-letter residue codes | `K143R`, `K128fs` |
| `REPORT` | formatted report line (see below) | `ERG11: "lanosterol 14-alpha demethylase" (missense_variant c.428A>G p.Lys143Arg; T:0 C:562)` |

The header line is `#`-prefixed (`#GENE  HGVSc  HGVSp  NT  AA  REPORT`), and a
`NA` denotes fields that have no HGVS string. The TSV prints to stdout unless
`--output` is given.

Each report line has the form `label: "product" (consequence HGVS.c HGVS.p; depths)`:

- **label**: the `--query_genes` term that matched the feature, or the `--bedfile`
  name if there is no query term.
  If there is no query list or match, the label is the
  normalized product name
  (e.g. `lanosterol.14-alpha.demethylase: "lanosterol 14-alpha demethylase" (...)`).
- **product**: the CDS product from the reference GFF. It is quoted so commas
  do not break the comma-delimited report.
- **consequence**: SnpEff's consequence term(s), `&`-joined, followed by the HGVS
  changes without their transcript and protein prefixes.
- **depths**: the read depth for each allele, read off the VCF's sample `AD`
  (else `RO`/`AO`) field; omitted when the record carries neither.

An entry is tied back to the annotation by its variant's coordinates, which
select every overlapping annotation unit — the transcripts where the annotation
has them, else genes, else bare CDS records. Its `Feature_ID` is consulted only to
choose between several units overlapping one variant (matched against the unit's
`ID`, `Name` or `transcript_id`), since SnpEff emits a separate entry per
transcript and each entry's HGVS strings belong to exactly one of them.

An entry is dropped when any of its consequence terms is suppressed, when it is
not annotated against a transcript (e.g. `intergenic_region`), when it carries
neither an HGVS.c nor an HGVS.p string, or when its variant resolves to no CDS
product — because nothing overlaps it, because several units do and none answers
to the entry's `Feature_ID`, or because the resolved unit carries no
`--feature_qualifier` attribute. A transcript SnpEff flags with
`WARNING_TRANSCRIPT_NO_START_CODON` or `WARNING_TRANSCRIPT_INCOMPLETE` is treated
as noncoding, so its protein change is reported as `NA`.

#### Deviations from HGVS

HGVS strings are SnpEff's, so the columns follow the
[HGVS recommendations](https://hgvs-nomenclature.org/) except where noted below.
One-letter amino acid codes, `*` for a stop codon, and the short frameshift form
(`p.Lys128fs`) are permitted by HGVS and are not deviations.

| deviation | HGVS | reported | applies to |
| --- | --- | --- | --- |
| reference sequence identifier dropped (the gene label stands in for it) | `NM_000001.1:c.428A>G` | `c.428A>G` | `NT`, `AA`, `REPORT` |
| parentheses around predicted protein changes absent | `p.(Lys143Arg)` | `p.Lys143Arg` | `HGVSp`, `AA`, `REPORT` |
| synonymous change repeats the reference residue instead of using `=` | `p.Asp164=` | `p.Asp164Asp`, `D164D` | `HGVSp`, `AA`, `REPORT` |
| deleted/duplicated bases listed | `c.383del` | `c.383delA` | `HGVSc`, `NT`, `REPORT` |
| protein-level duplication described as an insertion | `p.Gly294_Ser297dup` | `p.Ala292_Ser293insSerGlySerAla` | `HGVSp`, `AA`, `REPORT` |
| protein-level indel not shifted 3′ through a repeat | `p.Gln444_Gly448del` | `p.Phe432_Gly436del` | `HGVSp`, `AA`, `REPORT` |
| `c.`/`p.` coordinate prefix dropped | `c.428A>G`, `p.Lys143Arg` | `428A>G`, `K143R` | `NT`, `AA` |

```bash
theiagene report_variants \
  --vcf sample.snpeff.vcf \
  --reference_gff reference.gff

theiagene report_variants \
  --vcf sample.snpeff.vcf \
  --reference_gff reference.gff \
  --suppress synonymous_variant \
  --output VARIANT_REPORT.tsv
```

#### prepare_snpeff

Write the GFF and config that `snpeff build` needs to build a database for a
custom genome, `--genome_id`:

| output | content |
| --- | --- |
| `<data_dir>/<genome_id>/genes.gff` | the reference GFF's annotation section, without any embedded `##FASTA` section |
| `--output` (default `snpEff.config`) | a copy of `--template_config` (e.g. the config installed alongside `snpEff.jar`) with `data.dir` pointed at `--data_dir`, the genome registered under `--organism`, and its codon table(s) assigned |

SnpEff reads the genome's sequences from `<data_dir>/<genome_id>/sequences.fa`
rather than the GFF, so the reference FASTA is staged there separately.

Each contig is assigned the codon table its GFF features declare in their
`transl_table` attribute (an NCBI translation table number), so a nuclear genome
and its mitochondrion can each carry their own. Features declaring no
`transl_table` are ignored, and a contig none of whose features declare one
falls back to SnpEff's default (Standard) table — with a warning when no contig
declares one at all. `--translation_table` takes precedence over the GFF and
applies one table genome-wide; it is required when a contig's features declare
more than one table, which raises an error rather than guessing. An error, or a
table SnpEff does not support, leaves nothing written.

```bash
theiagene prepare_snpeff \
  --reference_gff reference.gff \
  --template_config /snpEff/snpEff.config \
  --data_dir snpeff_data \
  --genome_id cauris \
  --organism "Candidozyma auris"
zcat -f reference.fasta.gz > snpeff_data/cauris/sequences.fa

snpeff build -c snpEff.config -gff3 cauris
snpeff ann -c snpEff.config cauris sample.vcf > sample.snpeff.vcf
theiagene report_variants \
  --vcf sample.snpeff.vcf \
  --reference_gff snpeff_data/cauris/genes.gff
```

## Library

The subcommands share a gene/feature data model — the `Feature` and
`FeatureCol` classes — that turns a flat GFF3 annotation into a
navigable gene → RNA → CDS/exon hierarchy. See
[src/theiagene/lib/README.md](src/theiagene/lib/README.md) for a human-readable
introduction and full API reference.

Reading references is shared through `theiagene.lib.parsers`:

| function | returns |
| --- | --- |
| `assimilate_gff(gff)` | a `FeatureCol` of every GFF3 record, grouped into its hierarchy |
| `iter_gff_features(gff)` | each GFF3 record as a `Feature`, ungrouped |
| `iter_gff_lines(gff)` | the raw lines of a GFF3's annotation section |
| `gff_translation_tables(gff)` | a `{contig: table}` map of the `transl_table` each contig's features declare; raises `ValueError` for a contig declaring more than one |
| `import_vcf(vcf)` | a `pysam.VariantFile`, with GQ values written in scientific notation scrubbed to integers |
| `import_bam(bam)` | an indexed `pysam.AlignmentFile` |

A `.gz` GFF3 is read through `gzip`, and every GFF3 reader stops at an embedded
`##FASTA` directive (in any case, with or without a space after the `##`).