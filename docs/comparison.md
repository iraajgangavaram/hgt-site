# How hgt-site compares with other tools

Finding foreign DNA inserted into a host genome is an old problem with several
existing tools, and some of them are better than hgt-site at the thing they were
built for. This page is here so that someone can decide what to use, including
deciding not to use hgt-site.

Dates were read from the repositories on 2 November 2026.

| Tool | What it looks for | Latest release | Last commit |
| --- | --- | --- | --- |
| [LGTSeek](https://github.com/rileydavidr/lgtseek) | Bacterial lateral gene transfer in host sequencing data | v1.08, January 2015 | June 2016 |
| [IGS lgtseek-pipeline](https://github.com/IGS/lgtseek-pipeline) | The same method, packaged as a pipeline | none tagged | February 2022 |
| [VIRUSBreakend](https://github.com/PapenfussLab/gridss) | Viral integration, using single breakends | v2.13.2, February 2022 | May 2025 |
| [ViFi](https://github.com/namphuon/ViFi) | Viral integration, using HMMs over viral families | none tagged | July 2021 |
| [nf-core/hgtseq](https://github.com/nf-core/hgtseq) | A pipeline that aligns reads and classifies the unmapped ones | 1.1.0, May 2023 | actively developed |
| hgt-site | Bacterial or viral insertion sites, tested against chance chimeras | 0.4.0, November 2026 | actively developed |

## When another tool is the better choice

**You are looking for viral integration and want the best available evidence.**
Use [VIRUSBreakend](https://doi.org/10.1093/bioinformatics/btab343). It assembles
the breakend itself rather than inferring a position from read pairs, it handles
the viral taxonomy and the repetitive viral genomes that cause trouble, it is
part of GRIDSS and so benefits from a mature structural variant caller, and it
has been published and used on real cancer genomes. hgt-site has not been
compared with it. Until that comparison exists, anyone who tells you hgt-site is
better for viral integration, including this documentation, is guessing.

**You want horizontal gene transfer between bacterial genomes over evolutionary
time.** None of the tools above address that. It is a phylogenetic question and
needs phylogenetic methods.

**You have long reads.** A single read spanning a junction is stronger evidence
than anything a paired-read method can produce.

**You want a whole workflow rather than one step.** A pipeline will align your
reads, classify the unmapped ones, run quality control and assemble a report.
hgt-site does one step of that and expects the inputs to exist already.

## What hgt-site does that the others do not

**It tests each candidate against the number of chance chimeras expected
there.** The tools above filter on read support, on alignment evidence, or on
assembly of the junction. None of them model the rate at which library
preparation joins a host fragment to a foreign fragment by chance, which is the
dominant source of false sites. On the benchmark datasets, filtering on read
support alone reported 24,861 sites of which 390 were real.

**It works on any foreign taxon, not a curated reference set.** VIRUSBreakend
and ViFi both work against viral references. hgt-site takes whatever Kraken2
reported, so a bacterium nobody expected is treated the same as one that was.

**It keeps the rejected candidates.** Every candidate position stays in the
output with the reason it failed, so a region can be checked without rerunning
anything.

**It is small and installs in one line.** It is on Bioconda and PyPI, with a
container built from the Bioconda recipe. One dependency, no Java, no reference
bundle, no workflow engine. VIRUSBreakend needs GRIDSS and a Java runtime, and
ViFi needs its own reference data. For LGTSeek in particular, which has had no commit
since 2016, getting it to run at all now takes real work.

## What has actually been compared

Only LGTSeek, and only its read classification step rather than its later
BLAST-based filtering. That comparison is in the [benchmarks](benchmarks.md),
with the scripts and raw results in the repository.

Comparing hgt-site with VIRUSBreakend on viral integration data is the obvious
next benchmark and has not been done. It is a harder test than the LGTSeek
comparison and may well go against hgt-site.

## The honest summary

hgt-site is a focused tool with a specific statistical idea behind it, tested
carefully on simulated data and not at all on real data. VIRUSBreakend is a
mature, published tool for the viral case. If your problem is viral integration
in a real sample today, use VIRUSBreakend. If your problem is bacterial
sequence in a host sample, or you want to see the reasoning behind each call,
hgt-site is worth half an hour.
