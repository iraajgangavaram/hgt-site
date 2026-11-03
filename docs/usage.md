# Usage

## One command, from an alignment

```bash
hgt-site run \
  --bam sample.bam \
  --kraken-db /path/to/kraken2-db \
  --host-taxid 9606 \
  --threads 8 \
  -o sample.hgt-site.tsv \
  --bed sample.hgt-site.bed \
  --vcf sample.hgt-site.vcf
```

`run` pulls the unmapped reads out of the alignment, classifies them with Kraken2, and calls the sites. It needs `kraken2` on PATH:

```bash
conda install -c conda-forge -c bioconda kraken2
```

A Kraken2 database is also needed. The standard database is large, around 70 GB; a smaller one works if it holds the taxa you care about, and the pre-built databases at [benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2) include smaller options.

The unmapped reads and the Kraken2 output go to a temporary folder that is deleted afterwards. Pass `--workdir reads/` to keep them, which is worth doing if you want to look at the classifications yourself or run the calling step again with different thresholds. Extra Kraken2 options can be passed through with `--kraken-arg`, for example `--kraken-arg --confidence=0.1`.

## When the reads are already classified

Many pipelines run Kraken2 on the unmapped reads themselves. Use `call` with their output, which skips the classification step. Several Kraken2 files can be given at once, which suits a pipeline that writes one per read set:

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

## Options

| Option | What it does |
|---|---|
| `--bed sites.bed` | Also writes the sites that passed every filter as a BED file, for genome browsers |
| `--vcf sites.vcf` | Also writes every candidate site as a VCF file, with `PASS` or the failing flags in the FILTER column |
| `--annotation genes.gtf` | Reports the gene each site falls in, or the nearest one. GTF or GFF3, gzipped or not, as supplied by Ensembl, GENCODE or NCBI |
| `--max-gene-distance 100000` | How far from a gene a site can be and still have it named (default 100000) |
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
