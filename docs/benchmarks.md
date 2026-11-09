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

## Viral insertions

The benchmark above inserts bacterial DNA. Viral integration is the other case people look for, and the donor genome there is much smaller, so the whole benchmark was repeated with a virus in place of *Escherichia coli*.

| | Sites reported | Of those, real | Insertions found |
| --- | ---: | ---: | ---: |
| hgt-site, bacterial donor | 364 | 364 | 360 of 360 |
| hgt-site, viral donor | 361 | 357 | 357 of 360 |
| Read support alone, viral donor | 24,320 | 384 | 360 of 360 |

Eighteen runs for each donor, three seeds across six settings of coverage and chimera rate. With a viral donor, precision is 0.989 with a 95% confidence interval of 0.972 to 0.996, against 1.000 for the bacterial donor. Read support alone gives 0.016 either way, because the chimeras come from the contaminating bacteria rather than from the donor, so the problem is the same size in both cases.

The errors appear where there are fewest chimeras, not most. At 2,000 and 10,000 chimeric pairs the viral runs report one false site each in four of the twelve runs; at 30,000 they report none. The expected count at a position scales with how many chance chimeras the sample contains, so when chimeras are rare the expected count is small and a modest cluster of pairs clears the test more easily. A 30 kb donor gives less sequence over which to estimate that background than a 4.6 Mb one, which is why the effect shows up with a virus and not with a bacterium.

**The virus used is a stand-in.** It is SARS-CoV-2, chosen because it is the right order of magnitude in size and base composition, and it does not integrate into host genomes. The benchmark tests whether the statistics hold when the donor genome is small, not whether hgt-site finds real viral integrations. A benchmark on a virus that does integrate, such as HPV or HBV, and a comparison against VIRUSBreakend, are both still to be done.

To reproduce it, build a Kraken2 database containing the host, the four contaminants and the viral genome, then run the same script with the virus as the donor:

```bash
python3 benchmark/run_benchmark.py \
    --host host.fa --donor virus.fa \
    --contaminants klebsiella.fa salmonella.fa pseudomonas.fa acinetobacter.fa \
    --kraken-db k2db --nodes k2db/taxonomy/nodes.dmp \
    --coverage 10 --chimera-pairs 30000 --seed 1 \
    --outdir viral_run
```

### Samples that contain free donor as well as integrated donor

The chance chimeras above were built from the contaminating bacteria, so the foreign half never carried the virus that was really inserted. A real sample is rarely so tidy. A transduced cell line or a virus-positive tissue sample holds a lot of free donor alongside the integrated copy, and that free donor makes chimeras with host DNA during library preparation. Those chimeras carry the same taxon as the genuine insertion, so nothing about the classification separates them.

`--donor-chimera-fraction` sets the share of chance chimeras taken from the donor rather than the contaminants. The runs below use the viral donor at coverage 10, three seeds each.

| Donor share of chimeras | Sites | Real | False | Precision | Recall |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 61 | 60 | 1 | 0.984 | 60 of 60 |
| 0.1 | 61 | 60 | 1 | 0.984 | 60 of 60 |
| 0.3 | 60 | 60 | 0 | 1.000 | 60 of 60 |
| 0.6 | 59 | 59 | 0 | 1.000 | 59 of 60 |
| 0.6, 30,000 chimeras | 56 | 56 | 0 | 1.000 | 56 of 60 |

Precision holds. In the last row about 17,900 of the 30,000 chimeras carried the integrated virus, around 3,660 positions had enough read support to look like sites, and none of the 56 reported sites was false. Read support alone reported all 3,660.

What gives way is sensitivity. Chimeras made from the donor are anchored pairs, so they raise the rate at which that taxon is expected anywhere in the genome, and a cluster of donor reads becomes less surprising. The insertions that were missed all had few supporting pairs and were rejected as NOT_ENRICHED. That is the test being cautious rather than wrong, but it has a practical consequence: in a sample heavy with free donor, a weakly supported site will be dropped. More coverage, or reads that cross the junction, is the way to recover it. Loosening the threshold is not.

## Uneven host coverage

The expected number of chance chimeras at a position assumes foreign pairs fall evenly across the host genome. A chimera's host half is drawn from the host DNA in the library, so a region present at several copies contributes chimeras in proportion to its copy number, and the even assumption reads the excess as an insertion. This is the case in any tumour with copy-number gains, in an amplified or low-input library, and in RNA-seq.

`--amplify START:LENGTH:FACTOR` puts a region of the host genome into the simulated sample at several copies, with both the extra host reads and the extra chance chimeras that follow from it. The runs below use a 500 kb region at six copies, coverage 10 and 10,000 chimeric pairs.

| Background | Runs | Sites | Real | False | Precision | Insertions found |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Genome-wide rate (the default) | 3 | 69 | 59 | 10 | 0.855 | 59 of 60 |
| `--local-background` | 3 | 59 | 59 | 0 | 1.000 | 59 of 60 |
| Genome-wide rate, no amplified region | 3 | 59 | 59 | 0 | 1.000 | 59 of 60 |
| `--local-background`, no amplified region | 3 | 59 | 59 | 0 | 1.000 | 59 of 60 |

Every one of the 10 false sites was inside the amplified region, with 6 to 9 supporting pairs against an expected count between 0.17 and 0.63. Measuring coverage locally removed all 10 and lost no true site. The insertion missed in one seed was missed by both, so it is not a cost of the option. With no amplified region the two reported exactly the same sites on all three seeds.

The expectation is only ever raised, never lowered. A 10 kb bin at this coverage holds a few hundred reads, so the ratio carries a few per cent of noise, and an early version that also scaled downwards let a marginal candidate through on one of the even-coverage seeds. A region above the genome mean is therefore held to a higher bar, while a quiet region is tested against the genome-wide rate, which is the stricter of the two. The option can cost sensitivity and cannot cost precision.

This matters for the viral case in particular. An integrated virus often sits in a region that is itself amplified, and a tumour genome is rarely at an even copy number anywhere.

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


The simulated datasets, the known insertion sites and the full results for all 18 runs are packaged as a single archive, with one ready-to-run dataset included so another tool can be scored on the same data without rebuilding anything: [hgt-site-benchmark-v1.0.0.tar.gz](https://github.com/iraajgangavaram/hgt-site/releases/download/0.3.0/hgt-site-benchmark-v1.0.0.tar.gz). Its README explains the design and how to score a different tool against it.

## Limitations

- The benchmarks use real sequence but simulated reads, a Kraken2 database far smaller than a standard one, and chimeras placed at random positions. Real data may include chimera hotspots and more reads classified above species level.
- The LGTSeek comparison covers its read classification step only, not its later BLAST-based filtering.
- Sites are tested one taxon at a time, so an insertion whose reads Kraken2 splits between related species can lose support. Grouping Kraken2 calls at genus or family level before testing was tried and did not help: genus grouping changed nothing in the benchmark, and family grouping merged real insertions with contamination from related bacteria, which lowered recall.
- Fragment sizes are measured from properly paired reads, so if the input contains only pairs with an unmapped mate, the default clustering limits are used.
- Exact junctions need reads that cross the junction, so they are found for fewer sites at low coverage, and for none when the input holds only pairs with an unmapped mate.
- The speed comparison used five bacterial genomes as LGTSeek's reference, and the time for its second alignment depends on the size of that reference. All timings were measured on one 2-core machine and will differ on other hardware.
- Across every benchmark run, including earlier versions of the code, one false site appeared in 4 of 59 simulated datasets, all of them in the samples with the fewest chance chimeras. The raw results are in the repository.
- A single insertion is occasionally reported as two neighbouring sites. This happened to one insertion in 4 of 18 runs on fresh seeds.
