# hgt-site

hgt-site finds places where DNA from bacteria or other organisms may have been inserted into a host genome, such as a human genome, using sequencing data. Most bacterial reads in these samples come from contamination or from DNA fragments joined by chance in the lab, so hgt-site checks each possible site against what chance alone would produce and keeps only the sites with clearly more support.

It takes paired-end reads aligned to the host (BAM, CRAM or SAM) and Kraken2 classifications of the reads that did not align, and writes a table of candidate sites. It works with any host alignment and Kraken2 output, including the files produced by nf-core/hgtseq.

## Installation

```bash
pip install hgt-site
```

## Usage

```bash
hgt-site call \
  --bam sample.bam \
  --kraken sample.kraken2.txt \
  --host-taxid 9606 \
  --nodes taxonomy/nodes.dmp \
  -o sample.hgt-site.tsv
```

To try it on simulated data with known insertion sites:

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

## Documentation and benchmarks

The full documentation, the output format and the benchmark results, including a comparison with LGTSeek, are on GitHub: https://github.com/iraajgangavaram/hgt-site

## Citation

If you use hgt-site, please cite it: https://doi.org/10.5281/zenodo.23842712

## Licence

MIT
