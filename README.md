# hgt-site

[![tests](https://github.com/iraajgangavaram/hgt-site/actions/workflows/tests.yml/badge.svg)](https://github.com/iraajgangavaram/hgt-site/actions/workflows/tests.yml)

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
| `MIXED_TAXA` | Less than 80% of pairs share one taxon, and the site holds more pairs from other taxa than background chimeras would explain |
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
| `expected_other_pairs`, `other_pvalue` | Pairs from other taxa expected by chance, and the test for an excess of them |
| `flags`, `filter` | Flags raised and the overall filter result |

## Preliminary results on simulated data

The simulator builds a random host genome (two 100 kb chromosomes), inserts five 1.5 to 3 kb fragments of a donor genome, and generates 100 bp paired reads with known placement. It adds free contaminant read pairs and chimeric artefact pairs that join random host positions to contaminant sequence. Results over 30 random seeds each:

| Setting | Precision | Recall |
|---|---|---|
| 10x coverage, 2,000 contaminant pairs, 150 chimeras | 1.000 | 1.000 |
| 10x coverage, 4,000 contaminant pairs, 300 chimeras | 1.000 | 1.000 |
| 5x coverage, 2,000 contaminant pairs, 300 chimeras | 1.000 | 0.980 |
| 10x, heavy noise, clusters filtered on support only (no enrichment test) | 0.217 | 1.000 |

The last row shows why the enrichment test matters: filtering on read support alone accepts most chimera clusters.

## Benchmark on real genomes

A second benchmark uses real sequence. Twenty fragments (1 to 5 kb) of the *Escherichia coli* K-12 genome were inserted into a 10 Mb region of human chromosome 21 (chr21:20,000,000-30,000,000, from nf-core/test-datasets). Paired 150 bp reads were simulated with wgsim. The reads were mixed with 50,000 free contaminant read pairs from *Klebsiella pneumoniae*, *Salmonella enterica*, *Pseudomonas aeruginosa* and *Acinetobacter baumannii*, and with host-contaminant chimeric pairs of the kind formed during library preparation. Reads were aligned to the host with BWA-MEM and the unmapped reads were classified with Kraken2, using a database built from the host region and the five bacterial genomes.

hgt-site was compared with two approaches on the same data:

- Support only: the same clusters, accepted when they have at least 3 supporting pairs, with no enrichment test.
- LGTSeek: the read classification step of LGTSeek (`LGT::LGTSeek->bwaPostProcess` from [IGS/lgtseek-pipeline](https://github.com/IGS/lgtseek-pipeline)), with reads aligned to the host and to the five bacterial genomes. LGTSeek reports read pairs rather than sites, so its host-side reads were grouped with the same clustering settings and a site needed at least 3 pairs.

Results over 3 random seeds per setting (60 true insertions per row):

| Coverage | Chimeric pairs | hgt-site precision | hgt-site recall | LGTSeek precision | LGTSeek recall | Support-only precision | Support-only recall |
|---|---|---|---|---|---|---|---|
| 10x | 2,000 | 1.000 | 1.000 | 0.813 | 1.000 | 0.813 | 1.000 |
| 10x | 10,000 | 1.000 | 1.000 | 0.053 | 1.000 | 0.053 | 1.000 |
| 10x | 30,000 | 1.000 | 1.000 | 0.007 | 1.000 | 0.007 | 1.000 |
| 5x | 2,000 | 1.000 | 0.983 | 0.781 | 0.817 | 0.813 | 1.000 |
| 5x | 10,000 | 1.000 | 0.983 | 0.047 | 0.883 | 0.053 | 1.000 |
| 5x | 30,000 | 1.000 | 0.983 | 0.007 | 0.983 | 0.007 | 1.000 |

As chimeras become more common, approaches based on read support alone report thousands of candidate sites, of which fewer than 1% are real. hgt-site kept every call correct across all settings. The one insertion it missed at 5x coverage came from a region shared between *E. coli* and *K. pneumoniae*, where Kraken2 split the reads between the two species.

The scripts are in [`benchmark/`](benchmark/). `run_benchmark.py` runs one replicate and `summarise.py` builds the table from a folder of runs.

## Limitations

The simulator models alignment and Kraken2 classification directly rather than running an aligner and classifier, and its genomes are random sequence without repeats. The real-genome benchmark uses real sequence but simulated reads, a Kraken2 database far smaller than a standard one, and chimeras placed at random positions. Real data may include chimera hotspots and more reads classified above species level. The LGTSeek comparison covers its read classification step only, not its later BLAST-based filtering. Sites are tested one taxon at a time, so an insertion whose reads Kraken2 splits between related species can lose support. The default clustering settings assume fragments of around 350 bp.

## Future work

- Test on public sequencing data, including the datasets used to evaluate nf-core/hgtseq
- Group reads from closely related taxa before testing, so that split classifications do not reduce support
- Estimate fragment size from the data and set clustering limits from it
- Use split reads that cross the junction to refine breakpoints to base-pair resolution
- Add VCF and BED output
- Package for Bioconda and contribute an nf-core module

## Testing

```bash
pip install ".[test]"
pytest
```

## Use of generative AI

hgt-site is led by Iraaj Gangavaram, who chose the project from the options considered, set its scope, and made the decisions on its programming language, licence, benchmarking and comparison with LGTSeek. He reviewed and approved every change before it was committed. Claude (Anthropic), a generative AI assistant (Claude Opus 5.5), was used to write most of the source code, tests and benchmark scripts, to draft the documentation including this README, and to run and summarise the benchmarks. The methods were checked with the automated tests and the benchmarks reported above.

## Citation

If you use hgt-site in your work, please cite it. Use the "Cite this repository" button on the GitHub page, or the details in [CITATION.cff](CITATION.cff).

## Author

Iraaj Gangavaram ([ORCID 0009-0002-4615-2846](https://orcid.org/0009-0002-4615-2846))

## Licence

MIT. See [LICENSE](LICENSE).
