# hgt-site

hgt-site finds places where DNA from bacteria or other organisms may have been inserted into a host genome, such as a human genome, using sequencing data. Most bacterial reads in these samples come from contamination or from DNA fragments joined by chance in the lab, so hgt-site checks each possible site against what chance alone would produce and keeps only the sites with clearly more support.

It works with any host alignment and Kraken2 output, including the files that nf-core/hgtseq already produces.

<figure markdown>
  ![Precision against the number of chimeric read pairs](images/benchmark_precision_light.svg#only-light){ width="640" }
  ![Precision against the number of chimeric read pairs](images/benchmark_precision_dark.svg#only-dark){ width="640" }
  <figcaption>As chance chimeras become more common, read support alone reports thousands of false sites, while every hgt-site call stays a real insertion. Details in <a href="benchmarks/">Benchmarks</a>.</figcaption>
</figure>

<figure markdown>
  ![A terminal showing hgt-site being installed and run](images/demo_light.svg#only-light){ width="760" }
  ![A terminal showing hgt-site being installed and run](images/demo_dark.svg#only-dark){ width="760" }
</figure>

## In short

- Takes paired-end reads aligned to the host (BAM, CRAM or SAM) and Kraken2 calls for the reads that did not align.
- Groups the reads that point at the same host position, and tests each group against the number of chance chimeras expected there.
- Gives the junction to the base pair where reads cross it.
- Writes a table, a BED file, a VCF file and a MultiQC summary.
- Takes about 3 seconds on a sample of 2.1 million reads once the reads are aligned and classified.

## The problem

When a sample is aligned to its host genome, some reads do not map. Kraken2 often assigns many of these to bacteria or viruses. A small number may come from real integrations of foreign DNA into the host genome, but most come from contamination or from chimeric fragments formed during library preparation, where a host fragment and a contaminant fragment are joined by chance.

A pair where one mate maps to the host and the other is foreign therefore does not show an integration on its own. The useful signal is several such pairs, from the same foreign taxon, pointing at the same host position.

## Getting started

On your own data, with a host alignment and a Kraken2 database:

```bash
pip install hgt-site
hgt-site run --bam sample.bam --kraken-db /path/to/kraken2-db --host-taxid 9606 -o sites.tsv
```

Or try it in a few seconds on simulated data with known insertion sites, which needs nothing else installed:

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

[Installation](installation.md){ .md-button } [Usage](usage.md){ .md-button }
