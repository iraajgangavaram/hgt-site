# Benchmarks

Three things were measured: whether the sites hgt-site reports are real, how close its junctions are to the truth, and how fast it runs. The scripts are in [`benchmark/`](https://github.com/iraajgangavaram/hgt-site/tree/main/benchmark) and the numbers behind each figure are in `benchmark/results/`.

## Accuracy on real genomes

Twenty fragments (1 to 5 kb) of the *Escherichia coli* K-12 genome were inserted into a 10 Mb region of human chromosome 21 (chr21:20,000,000-30,000,000, from nf-core/test-datasets). Paired 150 bp reads were simulated with wgsim. The reads were mixed with 50,000 free contaminant read pairs from *Klebsiella pneumoniae*, *Salmonella enterica*, *Pseudomonas aeruginosa* and *Acinetobacter baumannii*, and with host-contaminant chimeric pairs of the kind formed during library preparation. Reads were aligned to the host with BWA-MEM and the unmapped reads were classified with Kraken2, using a database built from the host region and the five bacterial genomes.

hgt-site was compared with two approaches on the same data:

- **Support only**: the same clusters, accepted when they have at least 3 supporting pairs, with no enrichment test.
- **LGTSeek**: the read classification step of LGTSeek (`LGT::LGTSeek->bwaPostProcess` from [IGS/lgtseek-pipeline](https://github.com/IGS/lgtseek-pipeline)), with reads aligned to the host and to the five bacterial genomes. LGTSeek reports read pairs rather than sites, so its host-side reads were grouped with the same clustering settings and a site needed at least 3 pairs.

<figure markdown>
  ![Precision and recall against the number of chimeric read pairs](images/benchmark_accuracy_light.svg#only-light){ width="720" }
  ![Precision and recall against the number of chimeric read pairs](images/benchmark_accuracy_dark.svg#only-dark){ width="720" }
  <figcaption>Six random seeds per setting, three with each simulator, 20 true insertions in each. Small points are single runs and lines join the means. Precision is on a log scale, because the other two approaches fall below 1% while hgt-site stays at 100%.</figcaption>
</figure>

Results over 3 random seeds per setting, 60 true insertions per row:

| Coverage | Chimeric pairs | hgt-site precision | hgt-site recall | LGTSeek precision | LGTSeek recall | Support-only precision | Support-only recall |
|---|---|---|---|---|---|---|---|
| 10x | 2,000 | 1.000 | 1.000 | 0.813 | 1.000 | 0.813 | 1.000 |
| 10x | 10,000 | 1.000 | 1.000 | 0.053 | 1.000 | 0.053 | 1.000 |
| 10x | 30,000 | 1.000 | 1.000 | 0.007 | 1.000 | 0.007 | 1.000 |
| 5x | 2,000 | 1.000 | 0.983 | 0.781 | 0.817 | 0.813 | 1.000 |
| 5x | 10,000 | 1.000 | 0.983 | 0.047 | 0.883 | 0.053 | 1.000 |
| 5x | 30,000 | 1.000 | 0.983 | 0.007 | 0.983 | 0.007 | 1.000 |

As chimeras become more common, approaches based on read support alone report thousands of candidate sites, of which fewer than 1% are real. hgt-site kept every call correct across all settings. The one insertion it missed at 5x coverage came from a region shared between *E. coli* and *K. pneumoniae*, where Kraken2 split the reads between the two species.

### Independent simulator

To check that the results do not depend on how our own script places insertions, the benchmark was repeated with the integration simulator of Scott et al. (2022), from [aehrc/vector-integration-simulation-pipeline](https://github.com/aehrc/vector-integration-simulation-pipeline). Its insertion model adds rearranged and partly deleted fragments, deletions of host sequence at the site, and small gaps or overlaps at the junctions, and it places insertions with its own position model. It was run with *E. coli* fragments of 1 to 5 kb, no whole-genome insertions and no free episomes. Reads, contamination, chimeras, alignment and classification were the same as above, on 3 new random seeds per setting:

| Coverage | Chimeric pairs | hgt-site precision | hgt-site recall | LGTSeek precision | LGTSeek recall | Support-only precision | Support-only recall |
|---|---|---|---|---|---|---|---|
| 10x | 2,000 | 1.000 | 1.000 | 0.811 | 1.000 | 0.811 | 1.000 |
| 10x | 10,000 | 1.000 | 1.000 | 0.053 | 1.000 | 0.053 | 1.000 |
| 10x | 30,000 | 1.000 | 1.000 | 0.006 | 1.000 | 0.006 | 1.000 |
| 5x | 2,000 | 1.000 | 1.000 | 0.791 | 0.883 | 0.811 | 1.000 |
| 5x | 10,000 | 1.000 | 0.983 | 0.050 | 0.933 | 0.053 | 1.000 |
| 5x | 30,000 | 1.000 | 0.967 | 0.006 | 0.983 | 0.006 | 1.000 |

The pattern matches the first benchmark. The three insertions hgt-site missed at 5x coverage were again fragments whose reads Kraken2 split between *E. coli* and *K. pneumoniae*.

The simulator needs Biopython 1.77. For these runs, the few Biopython functions it uses (FASTA reading, sequence slicing and reverse complement) were supplied by a small compatibility module, which passed the simulator's own test suite.

## Libraries with long fragments

Reads from one insertion lie within about one fragment length of the junction, so with fixed clustering limits, libraries with long fragments can have one insertion split into two candidate sites. hgt-site measures the fragment size from properly paired reads and widens its limits when needed. The built-in simulator was run with mean fragment lengths of 400, 800 and 1000 bp, with 20 random seeds each, 5 insertions per seed, 300 chimeric pairs and 4,000 contaminant pairs:

<figure markdown>
  ![Calls per real insertion at three fragment lengths](images/benchmark_fragments_light.svg#only-light){ width="640" }
  ![Calls per real insertion at three fragment lengths](images/benchmark_fragments_dark.svg#only-dark){ width="640" }
</figure>

Both settings found every insertion. With fixed limits, most insertions in the 800 and 1000 bp libraries were reported twice, and precision at 1000 bp was 85%. With limits set from the measured fragment size, precision was 100% at 800 bp and 97.5% at 1000 bp. At 400 bp the limits are not changed and the results are identical.

## Junction accuracy

Read pairs place a junction only to within about one fragment length, so hgt-site also looks for reads that cross the junction and are soft-clipped past it. This was measured on 8 runs of the benchmark above (10x and 5x coverage, 2,000 and 30,000 chimeric pairs, 2 seeds each), using the full sorted alignment rather than one filtered down to pairs with an unmapped mate. Across the 157 sites that matched a true insertion, split reads were found for 154.

<figure markdown>
  ![How close each junction is, with and without split reads](images/benchmark_junctions_light.svg#only-light){ width="640" }
  ![How close each junction is, with and without split reads](images/benchmark_junctions_dark.svg#only-dark){ width="640" }
</figure>

With split reads, 71% of junctions were exact and 98% were within 2 bp, against 41% exact and 57% within 2 bp from read pairs alone. The remaining 2% are the three sites with no split reads, where the estimate from read pairs is kept. The same sites passed the filters either way, and the extra step made no measurable difference to the run time.

## Speed

hgt-site only needs read pairs where at least one mate is unmapped, and it skips every other read before doing any work on it. It reads the host alignment once and needs no second alignment.

Measured on one sample of 2.1 million paired reads (1 million pairs from human chr21:20,000,000-30,000,000 and 50,000 contaminant pairs, aligned with BWA-MEM), on 2 CPU cores:

| Step after the host alignment | hgt-site | LGTSeek |
|---|---|---|
| Kraken2 on the unmapped reads | 12 s | not needed |
| Second alignment of all reads to the bacterial genomes | not needed | 129 s |
| Calling sites or classifying reads | 2.7 s | 25 s |
| Total | 15 s | 154 s |

<figure markdown>
  ![Seconds after the host alignment for hgt-site and LGTSeek](images/benchmark_speed_light.svg#only-light){ width="720" }
  ![Seconds after the host alignment for hgt-site and LGTSeek](images/benchmark_speed_dark.svg#only-dark){ width="720" }
</figure>

hgt-site times are the median of five runs. Where a pipeline has already run Kraken2, hgt-site adds about 3 seconds for a sample of this size.

Run time grows in proportion to the size of the alignment file. The same sample was subsampled to about 0.5 and 1 million records and concatenated to 4.2 and 8.4 million records, and each file was run three times with each version:

<figure markdown>
  ![Run time against the number of records in the BAM file](images/benchmark_scaling_light.svg#only-light){ width="640" }
  ![Run time against the number of records in the BAM file](images/benchmark_scaling_dark.svg#only-dark){ width="640" }
</figure>

At 8.4 million records the current version took a median of 9.8 s, against 25.5 s for version 0.1.1, or about 0.9 million records a second. At that rate, a 100-million-read exome takes roughly two minutes. Times vary by about 10% between runs of the same file.

## Reproducing the benchmarks

`run_benchmark.py` runs one replicate, `summarise.py` builds the tables from a folder of runs, `junctions.py` measures junction accuracy, and `figures.py` draws the figures. Add `--intvi-script` with the path to the independent simulator's `scripts/insert_virus.py` to use it instead of the built-in insertion step, and `--lgtseek-lib` with LGTSeek's `perl5` library folder to include LGTSeek.

## Limitations

- The benchmarks use real sequence but simulated reads, a Kraken2 database far smaller than a standard one, and chimeras placed at random positions. Real data may include chimera hotspots and more reads classified above species level.
- The LGTSeek comparison covers its read classification step only, not its later BLAST-based filtering.
- Sites are tested one taxon at a time, so an insertion whose reads Kraken2 splits between related species can lose support. Grouping Kraken2 calls at genus or family level before testing was tried and did not help: genus grouping changed nothing in the benchmark, and family grouping merged real insertions with contamination from related bacteria, which lowered recall.
- Fragment sizes are measured from properly paired reads, so if the input contains only pairs with an unmapped mate, the default clustering limits are used.
- Exact junctions need reads that cross the junction, so they are found for fewer sites at low coverage, and for none when the input holds only pairs with an unmapped mate.
- The speed comparison used five bacterial genomes as LGTSeek's reference, and the time for its second alignment depends on the size of that reference. All timings were measured on one 2-core machine and will differ on other hardware.
- Across every benchmark run, including earlier versions of the code, one false site appeared in 4 of 59 simulated datasets, all of them in the samples with the fewest chance chimeras. The raw results are in the repository.
- A single insertion is occasionally reported as two neighbouring sites. This happened to one insertion in 4 of 18 runs on fresh seeds.
