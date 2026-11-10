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
| `gc_foreign`, `gc_difference`, `kmer_ratio` | Base composition of the foreign reads, and how far it sits from the host's own. `kmer_ratio` is near 1 when the reads are indistinguishable from host sequence, which usually means the classifier got them wrong |
| `depth_ratio` | Host coverage at the site over the genome mean, used to scale `expected_pairs`. Always 1 unless `--local-background` was given and the local estimate was trustworthy |
| `expected_other_pairs`, `other_pvalue` | Pairs from other taxa expected by chance, and the test for an excess of them |
| `flags`, `filter` | Flags raised and the overall filter result |

A site passes when `filter` is `PASS`. Otherwise `filter` lists the reasons it failed.

## BED

With `--bed`, the sites that passed every filter are written as BED: one line per site, named by the main taxid (for example `taxid_562`), with the number of supporting pairs as the score, capped at 1000 as the format requires. This loads straight into IGV or a genome browser.

## VCF

With `--vcf`, every candidate site is written as a VCF 4.2 insertion record, so the calls can go through the same tools as other structural variants.

- `POS` is the last host base before the insertion. It comes from the split reads where there are any, otherwise from the anchors.
- `ALT` is the symbolic allele `<INS>` and `SVTYPE=INS`. hgt-site does not read the reference sequence, so `REF` is `N`.
- `QUAL` is the q-value on the phred scale, so `QUAL` of 13 is a q-value of 0.05 and larger numbers mean stronger enrichment. It is capped at 1000. Sites that failed a filter still carry their score, so a reader can rank the failures.
- `FILTER` is `PASS` or the failing flags.
- INFO holds `TAXID`, `NPAIRS`, `NFWD`, `NREV`, `QVAL`, `SITE_START`, `SITE_END`, `SPLIT_LEFT`, `SPLIT_RIGHT`, and the `ONE_SIDED` flag.
- A junction estimated from read pairs rather than found in split reads carries `IMPRECISE`, with `CIPOS` and `CIEND` giving how far the real junction could be from `POS`. The interval runs between the estimates from the two sides of the junction, or across the whole region the anchors cover when only one side has any. A junction found by split reads is exact and carries no interval.
- There is one sample column, named after the alignment file or after `--sample-name`, holding `GT`, `PR` and `SR`. `PR` is the supporting read pairs and `SR` the split reads crossing the junction. `GT` is always `./.`, because hgt-site counts the reads supporting a junction and does not work out how many copies of the host carry it.
- The header records the exact command and every threshold used, as `##hgt_siteCommand` and `##hgt_siteSetting` lines, so a file can be traced back to how it was made. Contig lengths come from the alignment header.

To keep only the sites that passed, or to filter on the score:

```bash
bcftools view -f PASS sample.hgt-site.vcf
bcftools view -i 'QUAL>50' sample.hgt-site.vcf
```

If the name given to `--vcf` ends in `.gz`, the file is block compressed and indexed with tabix, so it can be queried by region and used directly in a pipeline:

```bash
hgt-site call --bam sample.bam --kraken classified.txt --host-taxid 9606 \
    -o sample.hgt-site.tsv --vcf sample.hgt-site.vcf.gz
bcftools view sample.hgt-site.vcf.gz chr8:128000000-129000000
```

## Genes

With `--annotation genes.gtf`, three more columns are filled in: `gene`, `gene_id` and `gene_distance`. The gene reported is the one containing the position hgt-site called, with `gene_distance` 0, or the nearest gene within `--max-gene-distance` when the site falls between genes, with the distance in base pairs. A site with no gene within that distance leaves all three blank rather than naming something half a megabase away.

Any GTF or GFF3 works, gzipped or not, including the files Ensembl, GENCODE and NCBI supply beside their reference genomes. Only the lines whose third column is `gene` are read, so a full annotation with every transcript and exon costs nothing extra. In the VCF the same information appears as `GENE` and `GENEDIST`.

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
| `EXCLUDED` | Overlaps a region given to `--exclude-regions` |
| `ONE_SIDED` | Anchors from only one side of the junction. Information only, and does not fail a site, since a long insertion can leave only one junction within insert range. |
