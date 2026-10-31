# hgt-site

hgt-site calls candidate horizontal gene transfer (HGT) sites in a host genome from paired-end sequencing data. It takes reads aligned to the host and Kraken2 classifications of the reads that did not align, finds places where one mate maps to the host and the other is classified as a foreign organism, and tests whether each place has more support than chance chimeras would give.

It is designed to run on the outputs that nf-core/hgtseq already produces.

## The problem

When a sample is aligned to its host genome, some reads do not map. Kraken2 often assigns many of these to bacteria or viruses. A small number may come from real integrations of foreign DNA into the host genome, but most come from contamination or from chimeric fragments formed during library preparation, where a host fragment and a contaminant fragment are joined by chance.

A pair where one mate maps to the host and the other is foreign therefore does not show an integration on its own. The useful signal is several such pairs, from the same foreign taxon, pointing at the same host position.

## How it works

1. Collect read pairs. The host alignment is read once. Pairs where one mate maps with mapping quality of at least 20 and the other mate is unmapped are kept as anchors. Pairs where both mates are unmapped are counted as background.
2. Keep foreign anchors. The unmapped mate's Kraken2 call is checked. Unclassified reads, the host taxon and uninformative taxa (root, cellular organisms) are removed. With an NCBI `nodes.dmp` file, the whole host lineage and the host's ancestors are removed too.
3. Cluster. Anchors are grouped by the end that faces the junction (the 3' end of a forward read, the 5' end of a reverse read), so reads from both sides of one insertion fall in the same cluster. A limit on cluster width stops scattered reads chaining together.
4. Test for enrichment. For each cluster, the number of pairs from its main taxon is compared with a Poisson expectation based on how many anchored pairs that taxon has across the whole genome. Chance chimeras are spread across the genome, while real integrations pile up at one place. P-values are adjusted with Benjamini-Hochberg, counting every genome window that could have produced a cluster.
5. Flag. Each site gets a filter value of `PASS` or a list of reasons it failed.

| Flag | Meaning |
|---|---|
| `LOW_SUPPORT` | Fewer than 3 supporting pairs |
| `MIXED_TAXA` | Less than 80% of pairs share one taxon |
| `LOW_MAPQ` | Mean mapping quality of the anchors below 20 |
| `LOW_COMPLEXITY` | Unmapped mates have low 3-mer entropy (repeats or low-complexity sequence) |
| `NOT_ENRICHED` | Q-value above 0.05 |
| `ONE_SIDED` | Anchors from only one side of the junction (information only, does not fail a site) |

All thresholds can be changed on the command line.

## Installation

```bash
git clone https://github.com/iraajgangavaram/hgt-site.git
cd hgt-site
pip install .
```

pysam is installed as a dependency and is needed for BAM and CRAM input. Plain SAM files can be read without it.

## Usage

```bash
hgt-site call \
  --bam sample.bam \
  --kraken sample.kraken2.txt \
  --host-taxid 9606 \
  --nodes taxonomy/nodes.dmp \
  -o sample.hgt-site.tsv
```

`--kraken` is the per-read output from Kraken2 (`--output`), run on the unmapped reads. Read names ending in `/1` or `/2` are handled, as is output made with `--use-names`.

To try it on simulated data with known insertion sites:

```bash
hgt-site simulate -o demo
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o demo.tsv
```

## Output

A tab-separated table with one row per candidate site:

| Column | Description |
|---|---|
| `chrom`, `start`, `end` | Host region covered by the anchors (0-based) |
| `left_junction`, `right_junction` | Estimated junction from forward and reverse anchors |
| `n_pairs`, `n_forward`, `n_reverse` | Supporting pairs, split by anchor strand |
| `top_taxid`, `top_taxid_fraction`, `n_taxa` | Main foreign taxon and how consistent the calls are |
| `mean_mapq` | Mean mapping quality of the anchors |
| `mean_mate_complexity` | Mean 3-mer entropy of the unmapped mates (0 to 1) |
| `taxon_anchored_pairs` | Anchored pairs for this taxon across the genome |
| `taxon_unanchored_pairs` | Pairs for this taxon where neither mate mapped |
| `expected_pairs`, `pvalue`, `qvalue` | Enrichment test |
| `flags`, `filter` | Flags raised and the overall filter result |

## Preliminary results on simulated data

The simulator builds a random host genome (two 100 kb chromosomes), inserts five 1.5 to 3 kb fragments of a donor genome, and generates 100 bp paired reads with known placement. It adds free contaminant read pairs and chimeric artefact pairs that join random host positions to contaminant sequence. Results over 30 random seeds each:

| Setting | Precision | Recall |
|---|---|---|
| 10x coverage, 2,000 contaminant pairs, 150 chimeras | 1.000 | 1.000 |
| 10x coverage, 4,000 contaminant pairs, 300 chimeras | 1.000 | 0.980 |
| 5x coverage, 2,000 contaminant pairs, 300 chimeras | 1.000 | 0.813 |
| 10x, heavy noise, clusters filtered on support only (no enrichment test) | 0.217 | 1.000 |

The last row shows why the enrichment test matters: filtering on read support alone accepts most chimera clusters.

## Limitations

The simulator models alignment and Kraken2 classification directly rather than running an aligner and classifier, and its genomes are random sequence without repeats. Real data will include repetitive regions, misclassified reads and chimera hotspots that the current simulation does not contain. The default clustering settings assume fragments of around 350 bp.

## Future work

- Benchmark on reads simulated from real genomes and aligned with BWA-MEM, and compare with LGTSeek and with the current nf-core/hgtseq output
- Estimate fragment size from the data and set clustering limits from it
- Use split reads that cross the junction to refine breakpoints to base-pair resolution
- Add VCF and BED output
- Package for Bioconda and contribute an nf-core module

## Testing

```bash
pip install ".[test]"
pytest
```

## Citation

If you use hgt-site in your work, please cite it. Use the "Cite this repository" button on the GitHub page, or the details in [CITATION.cff](CITATION.cff).

## Author

Iraaj Gangavaram ([ORCID 0009-0002-4615-2846](https://orcid.org/0009-0002-4615-2846))

## Licence

MIT. See [LICENSE](LICENSE).
