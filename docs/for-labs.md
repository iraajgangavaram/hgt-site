# Is hgt-site worth trying?

This page is for someone deciding whether to spend an afternoon on hgt-site, or
whether to let a student spend one. It is deliberately blunt about what the tool
has and has not been shown to do.

## When it is worth trying

- You have paired-end reads aligned to a host genome, and Kraken2 or a similar
  classifier assigns some of the unmapped reads to bacteria or viruses, and you
  want to know which of those positions are worth following up.
- You have a list of candidate integration sites from read support and it is too
  long to check by hand.
- You need the junction to the base pair so a primer can be designed around it.
- You want the reasoning behind each call to be visible, rather than a score from
  a model you cannot inspect.

## When it is not

- You want to detect horizontal gene transfer between bacterial genomes, over
  evolutionary time. That is a different problem and hgt-site does not address
  it. Phylogenetic methods do.
- You have long reads. hgt-site works on paired short reads, and a long read
  spanning a junction is better evidence than anything it can offer.
- You need a validated clinical result. See the limitations below.
- Your alignment has been filtered down to only the pairs with an unmapped mate.
  It will still run, but it cannot measure fragment size or find split reads at
  the junction, so the junctions stay approximate.

## What the evidence actually is

| Claim | What it rests on |
| --- | --- |
| Separates real insertions from chance chimeras | 18 simulated datasets, three seeds across six settings of coverage and chimera rate, built from real human and bacterial genomes. 364 sites reported, none false, all 360 insertions found. 95% confidence interval on precision 0.990 to 1.000. |
| Works for viral insertions as well as bacterial | Eighteen simulated datasets with a viral donor. 361 sites reported, 357 real, 357 of 360 insertions found, a precision of 0.989 against 1.000 for the bacterial donor. The virus used does not integrate and was chosen for its size, so this tests the statistics rather than the biology. |
| Better than filtering on read support | On the same datasets, read support alone reported 24,861 sites of which 390 were real. |
| Better than an existing tool | Benchmarked against LGTSeek, which reported thousands of false sites in the noisiest datasets. LGTSeek has had no commit since 2016. |
| Not an artefact of our own simulator | The benchmark was repeated with an independent published simulator and the conclusion held. |
| Junctions to the base pair | Over 8 runs and 157 insertions, 71% of junctions were exact and 98% within 2 bp, against 41% and 57% from read pairs alone. |
| Fast enough not to matter | About 3 seconds of calling on a sample of 2.1 million reads, on 2 cores. |
| Installs without a fight | On Bioconda and PyPI, with a container built from the Bioconda recipe. Bioconda builds and tests a package on its own infrastructure before publishing it, so the install has been checked on a clean machine by somebody other than the author. |
| Reproducible | Every benchmark rebuilds from a seed. The scripts and raw results are in the repository, and the datasets are archived separately. |

## What the evidence is not

**It has never been run on real sequencing data with a confirmed integration.**
Everything above is simulated reads from real genomes. Simulated reads are
cleaner than real ones in ways that matter here: the Kraken2 database used in the
benchmark holds six genomes rather than tens of thousands, chimeras are placed at
random rather than in hotspots, insertions are clean with no rearrangement at the
junction, and the host region is 10 Mb of one chromosome rather than a whole
genome with its repeats and segmental duplications.

A tool that does well on that benchmark has shown it can do the statistics. It
has not shown it works on a patient sample. Treat a site it reports as a lead to
follow up by PCR or by eye in a genome browser, not as a finding.

There is also no publication behind it. It is archived with a DOI and can be
cited, and the method is written out in full, but nobody outside the project has
reviewed it.

## A half-hour test

Run it on something where you already know the answer, before running it on
something where you do not.

```bash
conda install -c conda-forge -c bioconda hgt-site     # or: pip install hgt-site
hgt-site simulate -o demo --seed 1
hgt-site call --bam demo.sam --kraken demo.kraken2.txt --host-taxid 9606 -o sites.tsv
```

`demo.truth.bed` holds the insertions that were put in, so every site in
`sites.tsv` can be checked. There is also a
[notebook](https://colab.research.google.com/github/iraajgangavaram/hgt-site/blob/main/examples/quickstart.ipynb)
that does this in the browser.

Then run it on your own sample, and look at the rows that failed as well as the
rows that passed. Every candidate stays in the table with the reason it was
rejected, which is usually more informative than the sites that got through.

## If you do try it

A report of what happened is genuinely wanted, including a run that found
nothing. There is an
[issue template](https://github.com/iraajgangavaram/hgt-site/issues/new?template=real_data_report.yml)
for it that asks for the sample, the command and what came out, and nothing
about your data itself. Reports that go into the documentation are credited by
name.

Real data is the one thing this tool is missing, and it cannot be fixed from
inside the project.
