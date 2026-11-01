# Changelog

## Unreleased

- More than twice as fast on large BAM and CRAM files: reads whose pair is fully mapped are now skipped before any other work. A new `--threads` option uses more cores to read BAM and CRAM files.
- `--kraken` accepts several files, such as the two that nf-core/hgtseq writes for each sample.
- New `--bed` option writes the sites that passed every filter as a BED file.
- New `--vcf` option writes every candidate site as a VCF 4.2 insertion record, with PASS or the failing flags in FILTER.
- Clustering limits are widened automatically for libraries with fragments longer than about 500 bp, using fragment sizes measured from properly paired reads.
- hgt-site now reports how many anchored pairs it found and how many have a Kraken2 call, and warns if the Kraken2 file does not seem to match the alignment.
- `hgt-site simulate` now fills in mate positions and marks properly paired reads, so its SAM files convert cleanly to BAM.
- Python 3.10 or newer is now needed. Python 3.9 no longer gets security updates.

## 0.1.1 (2026-10-31)

- Published on PyPI. No changes to how hgt-site works.

## 0.1.0 (2026-10-31)

- First release.
