# Usage

```bash
hgt-site call \
  --bam sample.bam \
  --kraken sample.kraken2.txt \
  --host-taxid 9606 \
  --nodes taxonomy/nodes.dmp \
  -o sample.hgt-site.tsv
```

## Inputs

| Input | Notes |
|---|---|
| `--bam` | Paired-end reads aligned to the host, as BAM, CRAM or SAM. An indexed BAM or CRAM also lets hgt-site find exact junctions from split reads. |
| `--kraken` | The per-read output of Kraken2 (`--output`, not the report), run on the reads that did not map to the host. Read names ending in `/1` or `/2` are handled, as is output made with `--use-names`. Several files can be given at once. |
| `--host-taxid` | NCBI taxonomy ID of the host, for example 9606 for human or 10090 for mouse. |
| `--nodes` | Optional NCBI `nodes.dmp`. With it, reads classified anywhere in the host lineage are ignored, not only those assigned to the host taxid itself. |

## With nf-core/hgtseq outputs

hgtseq writes a sorted BAM and two Kraken2 per-read files for each sample, which is everything hgt-site needs:

```bash
hgt-site call \
  --bam results/preprocess/alignment/sample01_sorted.bam \
  --kraken results/results/classified/single_unmapped/sample01.kraken2.classifiedreads.txt \
           results/results/classified/both_unmapped/sample01.kraken2.classifiedreads.txt \
  --host-taxid 9606 \
  -o sample01.hgt-site.tsv \
  --bed sample01.hgt-site.bed \
  --vcf sample01.hgt-site.vcf \
  --multiqc sample01_mqc.json
```

## Options

| Option | What it does |
|---|---|
| `--bed sites.bed` | Also writes the sites that passed every filter as a BED file, for genome browsers |
| `--vcf sites.vcf` | Also writes every candidate site as a VCF file, with `PASS` or the failing flags in the FILTER column |
| `--multiqc sample_mqc.json` | Also writes a summary table that MultiQC picks up, so the counts appear in a pipeline report |
| `--sample-name` | Sample name for the MultiQC table (the alignment file name by default) |
| `--threads 4` | Uses more cores to read large BAM or CRAM files |
| `--reference genome.fa` | Reference FASTA, needed for CRAM input |
| `--exclude-taxid 1234` | Ignores a taxon, for example a known contaminant (can be repeated) |
| `--no-split-reads` | Skips the split-read search, so junctions come from read pairs only |

### Thresholds

Every threshold can be changed. `hgt-site call --help` lists them with their defaults.

| Option | Default | What it controls |
|---|---|---|
| `--min-mapq` | 20 | Minimum mapping quality for a read to count as an anchor |
| `--min-support` | 3 | Minimum supporting pairs for a site to pass |
| `--max-qvalue` | 0.05 | Maximum q-value of the enrichment test |
| `--min-taxon-fraction` | 0.8 | Minimum share of pairs from the main taxon before `MIXED_TAXA` is considered |
| `--max-other-pvalue` | 0.01 | P-value below which pairs from other taxa count as an excess |
| `--min-mean-mapq` | 20 | Minimum mean mapping quality of the anchors |
| `--min-complexity` | 0.5 | Minimum mean 3-mer entropy of the unmapped mates |
| `--max-gap`, `--max-span`, `--window` | 300, 800, 600 | Clustering limits in bp. Set any of them and the automatic widening for long fragments is turned off. |

## Messages while it runs

hgt-site reports how many anchored pairs it found, how many have a Kraken2 call, the fragment size it measured, and how many sites got an exact junction. If none of the anchored reads has a Kraken2 call, it warns that the Kraken2 file may not match the alignment. If the alignment has no index, it says that junctions come from read pairs only.

## Trying it without real data

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

This writes `demo.sam`, `demo.kraken2.txt`, `demo.truth.bed` and `demo.host.fa`. The BED file holds the five insertion sites the simulator used, so the calls can be checked against them.
