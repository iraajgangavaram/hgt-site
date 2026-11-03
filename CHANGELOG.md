# Changelog

## Unreleased

- New `--annotation` option takes a GTF or GFF3 file and reports the gene each site falls in, or the nearest gene within `--max-gene-distance`, in three new columns and in the VCF as `GENE` and `GENEDIST`. A position on its own is hard to act on; the gene is the first thing anyone asks for.
- The VCF now carries a `QUAL` score, the q-value on the phred scale, so the file can be filtered with `bcftools view -i 'QUAL>50'` instead of being parsed by hand.
- A junction estimated from read pairs now carries `CIPOS` and `CIEND`, saying how far the real junction could be from the reported position. A junction found in split reads is exact and carries no interval. This is the standard way to state positional uncertainty, and the SV comparison tools use it to decide whether two calls match.
- The VCF header records the exact command and every threshold used, as `##hgt_siteCommand` and `##hgt_siteSetting` lines, and `##contig` lines now carry contig lengths from the alignment header.
- The VCF has one sample column, named after the alignment or after `--sample-name`, with the supporting pair count in `PR` and the split-read count in `SR`. `GT` is always `./.`, because hgt-site counts reads supporting a junction and does not genotype it.
- A `--vcf` name ending in `.gz` is block compressed and indexed with tabix, so it can be queried by region.
- Fixed: the command recorded in the VCF header came from `sys.argv`, so calling `main()` from another program recorded that program's arguments.

## 0.3.0 (2026-11-02)

- New `hgt-site run` command: give it a host alignment and a Kraken2 database and it pulls out the unmapped reads, classifies them with Kraken2 and calls the sites in one step. `hgt-site call` is unchanged, for reads a pipeline has already classified.
- Junctions are now exact where reads cross them. Reads that are soft-clipped at the junction place the insertion point to the base pair, reported in four new columns and used for the VCF position. In the benchmark, 71% of junctions were exact and 98% were within 2 bp, against 41% and 57% from read pairs alone. This needs an indexed BAM or CRAM; without an index the estimate from read pairs is kept. `--no-split-reads` turns it off.
- New `--multiqc` option writes a summary table that MultiQC picks up, so the counts appear in a pipeline report.
- `hgt-site simulate` also writes the host reference, so its reads can be stored as CRAM.
- Every option now has help text of its own.
- Documentation site at https://iraajgangavaram.github.io/hgt-site/

## 0.2.0 (2026-11-01)

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
