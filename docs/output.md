# Output

## The table

One row per candidate site, tab separated.

| Column | Description |
|---|---|
| `chrom`, `start`, `end` | Host region covered by the anchors (0-based) |
| `left_junction`, `right_junction` | Junction estimated from forward and reverse anchors |
| `left_breakpoint`, `right_breakpoint` | Exact junction from split reads, where there are any |
| `left_split_reads`, `right_split_reads` | Reads clipped at that exact base |
| `n_pairs`, `n_forward`, `n_reverse` | Supporting pairs, split by anchor strand |
| `top_taxid`, `top_taxid_fraction`, `n_taxa` | Main foreign taxon and how consistent the calls are |
| `mean_mapq` | Mean mapping quality of the anchors |
| `mean_mate_complexity` | Mean 3-mer entropy of the unmapped mates (0 to 1) |
| `taxon_anchored_pairs` | Anchored pairs for this taxon across the genome |
| `taxon_unanchored_pairs` | Pairs for this taxon where neither mate mapped |
| `expected_pairs`, `pvalue`, `qvalue` | Enrichment test |
| `expected_other_pairs`, `other_pvalue` | Pairs from other taxa expected by chance, and the test for an excess of them |
| `flags`, `filter` | Flags raised and the overall filter result |

A site passes when `filter` is `PASS`. Otherwise `filter` lists the reasons it failed.

## BED

With `--bed`, the sites that passed every filter are written as BED: one line per site, named by the main taxid (for example `taxid_562`), with the number of supporting pairs as the score, capped at 1000 as the format requires. This loads straight into IGV or a genome browser.

## VCF

With `--vcf`, every candidate site is written as a VCF 4.2 insertion record, so the calls can go through the same tools as other structural variants.

- `POS` is the last host base before the insertion. It comes from the split reads where there are any, otherwise from the anchors.
- `ALT` is the symbolic allele `<INS>` and `SVTYPE=INS`. hgt-site does not read the reference sequence, so `REF` is `N`.
- `FILTER` is `PASS` or the failing flags.
- INFO holds `TAXID`, `NPAIRS`, `NFWD`, `NREV`, `QVAL`, `SITE_START`, `SITE_END`, `SPLIT_LEFT`, `SPLIT_RIGHT`, the `IMPRECISE` flag when no split reads were found, and the `ONE_SIDED` flag.

To keep only the sites that passed:

```bash
bcftools view -f PASS sample.hgt-site.vcf
```

## MultiQC

With `--multiqc sample_mqc.json`, hgt-site writes a custom-content table that MultiQC finds on its own, as long as the file name ends in `_mqc.json`. It holds the number of candidate sites, the number that passed, how many have an exact junction, and the read counts behind them.

## Flags

| Flag | Meaning |
|---|---|
| `LOW_SUPPORT` | Fewer than 3 supporting pairs |
| `MIXED_TAXA` | Less than 80% of pairs share one taxon, and the site holds more pairs from other taxa than background chimeras would explain |
| `LOW_MAPQ` | Mean mapping quality of the anchors below 20 |
| `LOW_COMPLEXITY` | Unmapped mates have low 3-mer entropy (repeats or low-complexity sequence) |
| `NOT_ENRICHED` | Q-value above 0.05 |
| `ONE_SIDED` | Anchors from only one side of the junction. Information only, and does not fail a site, since a long insertion can leave only one junction within insert range. |
