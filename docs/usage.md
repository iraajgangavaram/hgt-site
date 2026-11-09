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

## When you already know what you are looking for

If the foreign DNA is a vector, a transgene or a particular viral genome, there is no reason to build a classifier database to find out where it went. Give hgt-site the sequence instead and it matches the unmapped reads against it directly:

```bash
hgt-site call \
  --bam sample.bam \
  --donor vector.fa \
  --host-taxid 9606 \
  -o sample.hgt-site.tsv
```

No Kraken2, no database, no taxonomy files. Everything after that step is the same, so the enrichment test, the flags, the split-read junctions and all the output formats behave exactly as they do with classifier output.

Matching is exact 31-mer matching in both orientations, and a read needs five matching k-mers to count. A 150 bp read holds 120 k-mers and one sequencing error spoils at most 31 of them, so the bar is low for a real donor read and high for anything else. It runs at roughly 50,000 reads a second, or about 20 seconds per million unmapped reads.

Donor reads are reported under taxid 32630, NCBI's taxon for a sequence that was built rather than sampled, which is the honest label for a vector. Use `--donor-taxid` to report them under the real taxon if your donor has one, for example 333760 for HPV16.

Two things to know. This is k-mer matching and not alignment, so a donor that differs from your reads throughout, by more than a few per cent, will not be found this way: align the reads and use a classifier for that case. And a donor sharing sequence with the host, a vector carrying a human promoter for instance, will match host reads too, so mask those parts of the FASTA or expect the shared region to produce sites.

## When host coverage is uneven

The test asks how many foreign pairs chance would put at a position, and by default it assumes they fall evenly across the genome. A chimera is a host fragment joined to a foreign fragment, so the chance its host half comes from a given region depends on how much host DNA from that region is in the library. Three copies of a region contribute three times the chimeras, and the even assumption will read the excess as an insertion.

That matters in a tumour with copy-number gains, in an amplified or low-input library, and in RNA-seq where coverage follows expression over several orders of magnitude.

```bash
hgt-site call \
  --bam sample.bam \
  --classifications sample.kraken2.txt \
  --host-taxid 9606 \
  --local-background \
  -o sample.hgt-site.tsv
```

Host coverage is counted in 10 kb bins during the scan the tool already does, and the expected count at each site is multiplied by the local rate over the genome mean. That ratio is reported in the `depth_ratio` column, where 1 means the site was tested against the genome-wide rate as usual.

A site keeps a ratio of 1, and so is tested exactly as it would be without the option, when the local estimate cannot be trusted: fewer than 200 mapped reads in the bins around it, or a ratio below 0.1 or above 10, which is more often a mapping artefact or a reference gap than real copy number.

**It needs the whole alignment.** If you pass a BAM already filtered down to the pairs with an unmapped mate, the host coverage has been thrown away. hgt-site notices, says so, and falls back to the genome-wide rate rather than giving you a confident wrong answer.

On the benchmark, with a 500 kb region present at six copies, the even assumption reported 69 sites across three seeds of which 10 were false, all of them inside the amplified region. With `--local-background` the same three runs reported 59 sites and none was false, with no true site lost. On a run with no amplified region the two gave the same sites, since the measured ratio is then close to 1 everywhere.

## Several samples at once

Call each sample as usual, then merge the tables:

```bash
hgt-site cohort sample*.hgt-site.tsv -o cohort.tsv --matrix matrix.tsv
```

Each position is reported once, with the samples it appeared in. A site in most of the samples is flagged `COMMON`, which usually means the protocol, the reference or the reagents rather than the biology, and a site in more than one is flagged `RECURRENT`. `TAXA_DIFFER` means samples disagreed about which organism is there, which points at a mismapping-prone region.

This reads the tables rather than the alignments, so samples called at different times can be combined without rerunning anything. `--names` sets the sample names, `--window` changes how close two positions must be to count as the same site (500 bp by default), `--common-fraction` changes the share of samples needed for `COMMON`, and `--include-failed` merges every candidate rather than only the sites that passed.

## Inputs

| Input | Notes |
|---|---|
| `--bam` | Paired-end reads aligned to the host, as BAM, CRAM or SAM. An indexed BAM or CRAM also lets hgt-site find exact junctions from split reads. |
| `--classifications` | Per-read taxonomic calls for the reads that did not map to the host, also accepted as `--kraken`. Kraken2 and KrakenUniq (`--output`, not the report), Kaiju, Centrifuge, DIAMOND `--outfmt 102` and BLAST tabular output all work, as does any two-column file of read name and taxid. Gzipped files are fine. Read names ending in `/1` or `/2` are handled, as is Kraken2 output made with `--use-names`. Several files can be given at once. |
| `--classifier` | Which classifier wrote those files, if the automatic detection gets it wrong. One of kraken2, kaiju, centrifuge, diamond, blast, tsv. |
| `--donor` | FASTA of the donor sequence, gzipped or not, instead of `--classifications`. The unmapped reads are matched against it directly, so no classifier and no database are needed. `--donor-taxid` sets the taxid they are reported under. |
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
| `--local-background` | Works out how many chance chimeras to expect at a site from host coverage near that site, rather than from one genome-wide rate. Worth using wherever coverage is uneven: a tumour with copy-number changes, an amplified or low-input library, RNA-seq. Needs the whole alignment, not one filtered down to the pairs with an unmapped mate |
| `--exclude-regions artefacts.bed` | Leaves out positions listed in a BED file. A candidate overlapping one is reported with `EXCLUDED` instead of passing, for the artefacts a protocol or a reference produces in every sample |
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
