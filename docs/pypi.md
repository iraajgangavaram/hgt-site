# hgt-site

hgt-site finds places where DNA from bacteria or other organisms may have been inserted into a host genome, such as a human genome, using sequencing data. Most bacterial reads in these samples come from contamination or from DNA fragments joined by chance in the lab, so hgt-site checks each possible site against what chance alone would produce and keeps only the sites with clearly more support.

It takes paired-end reads aligned to the host (BAM, CRAM or SAM) and the per-read output of whatever classifier you ran on the reads that did not align, and writes a table of candidate sites, with optional BED, VCF and MultiQC files.

- **Any classifier.** Kraken2, KrakenUniq, Kaiju, Centrifuge, DIAMOND and BLAST tabular output are read directly, as is any two-column file of read name and taxid. All of them give identical calls on the same reads.
- **Junctions to the base pair** where reads cross the junction, rather than an estimate from read pairs.
- **The gene at each site**, with `--annotation` and a GTF or GFF3.
- **`hgt-site selftest`** proves the installation works on data with known answers, in a few seconds, with no database and no network.
- **One command** from a BAM and a Kraken2 database, with `hgt-site run`.

In benchmarks on real human and bacterial genomes with added contamination, hgt-site reported no false sites across 18 runs on fresh random seeds, 364 calls in all, while keeping every position with enough read support gave 24,861 sites of which 390 were real. On a sample of 2.1 million reads it took about 3 seconds once the reads were aligned and classified, against about 2.5 minutes for LGTSeek.

## Installation

```bash
pip install hgt-site
```

or, with conda, `conda install -c conda-forge -c bioconda hgt-site`.

## Usage

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

This classifies the unmapped reads with Kraken2 and calls the sites in one step, and needs `kraken2` on PATH. If the reads have already been classified, use `hgt-site call --kraken sample.kraken2.txt` instead, which accepts several files at once.

To try it on simulated data with known insertion sites:

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

## Documentation and benchmarks

Full documentation: https://iraajgangavaram.github.io/hgt-site/

Source code, benchmark scripts and results: https://github.com/iraajgangavaram/hgt-site

## Citation

If you use hgt-site, please cite it: https://doi.org/10.5281/zenodo.23842712

## Licence

MIT
