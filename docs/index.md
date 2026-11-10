# Which foreign reads in your host sample actually mean something?

You aligned your sample to its host genome, ran Kraken2 on the reads that did not map, and got back a pile of bacteria and viruses. Most of it is free contamination, or chimeric fragments made during library preparation, where a host fragment and a foreign fragment are joined at an arbitrary position. Those chimeras look exactly like integration sites and there are usually thousands of them.

hgt-site groups the read pairs that point at the same host position and tests each group against the number of pairs chance alone would put there, given how common that taxon is across the sample. What comes back is the positions with more support than contamination and chimeras can explain, and every rejected candidate with the reason it failed. Where reads cross the junction itself, the insertion point is given to the base pair.

It works on bacterial and viral taxa alike, and on a vector or transgene you supply the sequence of. It takes any host alignment with per-read calls from Kraken2, KrakenUniq, Kaiju, Centrifuge, DIAMOND or BLAST, and runs in seconds.

If the foreign DNA is something you put there, `--donor` takes its sequence directly and no classifier or database is needed at all.

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

- Takes paired-end reads aligned to the host (BAM, CRAM or SAM) and per-read taxonomic calls for the reads that did not align, from any of six classifiers, or the donor's own sequence instead.
- Groups the reads that point at the same host position, and tests each group against the number of chance chimeras expected there.
- Gives the junction to the base pair where reads cross it.
- Writes a table, a BED file, a VCF file and a MultiQC summary.
- Takes about 3 seconds on a sample of 2.1 million reads once the reads are aligned and classified.
- Handles bacterial and viral insertions, and vectors or transgenes given the construct's sequence.

## The problem

When a sample is aligned to its host genome, some reads do not map. Kraken2 often assigns many of these to bacteria or viruses. A small number may come from real integrations of foreign DNA into the host genome, but most come from contamination or from chimeric fragments formed during library preparation, where a host fragment and a contaminant fragment are joined by chance.

A pair where one mate maps to the host and the other is foreign therefore does not show an integration on its own. The useful signal is several such pairs, from the same foreign taxon, pointing at the same host position.

## Getting started

On your own data, with a host alignment and a Kraken2 database:

```bash
conda install -c conda-forge -c bioconda hgt-site     # or: pip install hgt-site
hgt-site run --bam sample.bam --kraken-db /path/to/kraken2-db --host-taxid 9606 -o sites.tsv
```

Or try it in a few seconds on simulated data with known insertion sites, which needs nothing else installed:

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

[Installation](installation.md){ .md-button } [Usage](usage.md){ .md-button }
