# Organism test data

Size-minimized subsets of three production `TheiaEuk_Illumina_PE` samples,
one per organism `medea_magic` runs gene coverage and SnpEff annotation for.
They are consumed by [`tests/test_organisms.py`](../test_organisms.py).

| directory | sample | query genes | codon table | variant caller |
| --- | --- | --- | --- | --- |
| `candidozyma_auris/` | 22-00052321 (Clade III) | `medea_magic` Clade III set (9) | 12 (CTG clade) | GATK |
| `aspergillus_fumigatus/` | SRR17660939 | cyp51A, AFUA_6G05300, AFUA_4G08340 | 1 | GATK |
| `cryptococcus_neoformans/` | SRR16002855 | `sterol.14-demethylase` (ERG11, CNA00300) | 1 | FreeBayes (Snippy) |

`samples.json` records each sample's query genes, the GCS sources and, for the
two samples run on the current `medea_magic`, the gene coverage Terra reported.

## Files

Each organism directory holds:

| file | contents |
| --- | --- |
| `reference.gff` | every GFF feature overlapping a query gene (expanded to whole overlapping genes), plus the contigs' `region` lines |
| `reference.fasta` | those gene spans' sequence, one record per span named `<contig>:<start>-<end>` (1-based) |
| `reads.bam(.bai)` | the alignments overlapping those spans; aux tags dropped and read names shortened |
| `variants.vcf` | the VCF records within 1 kb of those spans |
| `snpeff.vcf` | SnpEff 5.4c annotation of the full-data `extract_variants` output |
| `expected/coverage_stats.tsv` | `gene_coverage` over the full BAM and GFF |
| `expected/variant_report.tsv` | `report_variants` over the full-data annotation |
| `snpeff_unencoded.vcf` | *C. auris* only: the production annotation, made before `prepare_snpeff` percent-encoded the pipes in NCBI's `rna-gnl\|WGS:PGLS\|...` IDs |

## Provenance

- Samples come from the `theiaeuk_illumina_pe` table of
  `theiagen-training-workspaces/Theiagen_Konkel_Sandbox`.
- The full-data outputs were generated with this checkout's theiagene and
  SnpEff 5.4c. They follow the order and arguments of PHB
  `task_gene_coverage`/`task_snpeff`, and the `medea_magic` coverage thresholds
  (depth 10, mapping quality 60, base quality 13).
- The *A. fumigatus* report and the *A. fumigatus* and *C. auris* coverage are
  identical to the production Terra outputs. A samtools depth calculation over
  the query CDS gave the same means for all three samples.
- *C. neoformans* SRR16002855 predates the `medea_magic` GATK track. Its BAM/VCF
  are the Snippy outputs against the same GCF_000091045.1 assembly, with
  contigs renamed to the reference's versioned accessions (`NC_006670` ->
  `NC_006670.1`) and Snippy's own `ANN` removed before re-annotation.

## Regenerating

[`build_test_data.py`](build_test_data.py) rebuilds everything from the GCS
sources. It needs `gsutil` with workspace read access, pysam, samtools and
docker with a SnpEff 5.4c image. It fails unless the minimized data reproduce
the full-data outputs exactly.

```bash
python3 tests/data/build_test_data.py --work /path/to/scratch \
  --samtools "$(command -v samtools)" \
  --snpeff_image us-docker.pkg.dev/general-theiagen/theiagen/snpeff:5.4c
```

Regenerate the `expected/` outputs whenever a change to theiagene changes them
on purpose.
