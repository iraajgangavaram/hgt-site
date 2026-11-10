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
| [VISPA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/) | Vector integration sites from amplicon data, for gene therapy | v2, 2017 | see publication |
| hgt-site | Bacterial, viral or vector insertion sites, tested against chance chimeras | 0.4.0, November 2026 | actively developed |

## Viral integration

hgt-site and [VIRUSBreakend](https://doi.org/10.1093/bioinformatics/btab343) both
find viral DNA integrated into a host genome, and they go about it differently.

VIRUSBreakend assembles the junction itself with GRIDSS, builds a reference from
the viruses Kraken2 finds in the sample, and has been published and used on real
cancer genomes. It needs Java, the GRIDSS stack and its own virus database, and
runs in hours.

hgt-site tests each candidate position against the number of chance chimeras
expected there, takes the classifier output you already have, installs in one
line and runs in seconds. On simulated viral insertions it reports a precision
of 0.989.

The two have not been run on the same data, so there is no measured comparison
between them. For the strongest evidence on a single integration, assembling the
breakend is the more direct approach. hgt-site handles viral and bacterial
insertions through the same command and the same outputs, working from whatever
your classifier reported, which matters when you do not know in advance what the
foreign DNA is.

## Vector and transgene integration

Finding where a vector or a transgene landed is the same problem with a known donor, and it has its own tools. [VISPA2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5702242/) and similar packages are built for integration site analysis in gene therapy, working from amplicon data produced by LAM-PCR or sonication-linker PCR, where the library is deliberately enriched around the junctions.

That difference in data type matters more than any difference in method. Those tools expect enrichment and are built to count and track clones from it. hgt-site expects an ordinary whole-genome library and tests each candidate against the chimera rate in it, so an enriched library breaks its assumption outright: if the reads around junctions have been amplified on purpose, a genome-wide rate says nothing about what to expect at a junction.

So the two are not alternatives for the same data. Use an amplicon-based tool on amplicon data. hgt-site suits the case where you have plain whole-genome sequencing of a transduced or transgenic sample and want the insertion sites out of it, with `--donor` pointing at the construct and no classifier involved.

What has been measured is in the [benchmarks](benchmarks.md#vectors-and-transgenes), and it is simulated data only. No real transduced sample has been tried, and a real plasmid carries repeated elements and often a promoter shared with the host, neither of which the simulated construct had.

## When another tool is the better choice

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
