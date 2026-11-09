# Is it contamination or a real insertion?

You aligned a host sample, pulled out the reads that did not map, ran Kraken2
or something like it, and got bacterial or viral hits. Some of them sit at a
consistent position in the host genome. The question is whether that means
foreign DNA is inserted there, or whether it means nothing at all.

Usually it means nothing at all. This page is about how to tell, and most of
it applies whatever tool you use.

## Where foreign reads come from

**Contamination in the sample or the reagents.** Kits, buffers and surfaces
carry bacterial DNA, and low-biomass samples are mostly reagent in places.
These reads have no fixed relationship to the host genome, so they turn up
everywhere.

**Chimeras made during library preparation.** A host fragment and a foreign
fragment get joined into one molecule, by ligation or during amplification.
Sequenced, that molecule looks exactly like an insertion junction: one mate on
the host, one mate on something foreign, consistent with each other. There is
nothing about a single read pair that distinguishes a chimera from the real
thing.

**Misclassification.** A host read assigned to a bacterium because it falls in
a conserved or repetitive region. Common with small Kraken2 databases and with
anything repeat-rich.

**An actual insertion.** Real, and rarer than the first three by a long way.

## What does not settle it

**Read count on its own.** Pick a support threshold, say three pairs, and a
whole-genome sample will hand you thousands of positions. In a simulated
benchmark on 10 Mb of human sequence, every position with at least three
supporting pairs gave 24,861 candidates of which 390 were real. Raising the
threshold loses true sites faster than it loses false ones, because chimeras
accumulate wherever coverage is high.

**A convincing pileup in a browser.** Chimeras cluster too. Coverage spikes,
mapping hotspots and segmental duplications all produce tidy-looking stacks of
foreign-classified reads.

**The taxon being biologically plausible.** If you are looking at a gut sample,
every chimera in it will also involve a plausible gut organism.

## What does

**How much support there is against how much you would expect.** This is the
one that matters. Count the foreign reads from that taxon across the whole
genome, work out the rate per base, and ask how many you would expect in a
window the size of your candidate. A cluster of six pairs is remarkable in a
clean sample and unremarkable in one riddled with chimeras. Without this
comparison you cannot interpret a read count at all.

**Whether the support comes from one organism.** A real insertion is one
donor. A mixture of genera at one position is the signature of a chimera
hotspot, because chance does not care which contaminant it joins to the host.

**Reads that cross the junction.** A read soft-clipped at the breakpoint, with
the clipped part matching the donor, places the junction to the base pair and
is far stronger evidence than any number of read pairs. If you have them, they
settle it.

**Whether the position appears in every sample.** A site in all of your
samples, including the ones that should be negative, is almost always the
protocol or the reference rather than biology.

## Doing this with hgt-site

hgt-site does the four checks above on the files you already have. It takes the
host alignment and the per-read classifications, clusters the anchored pairs
into candidate positions, and tests each one against the number of chance
chimeras expected there.

```bash
conda install -c conda-forge -c bioconda hgt-site

hgt-site call \
    --bam host.bam \
    --classifications sample.kraken2.txt \
    --host-taxid 9606 \
    -o sites.tsv
```

Classifier output from KrakenUniq, Kaiju, Centrifuge, DIAMOND and BLAST works
as well, and so does a plain two-column file of read name and taxid. The format
is worked out from the file.

Every candidate stays in the table with the reason it failed, which is more
informative than the list of sites that passed:

- `NOT_ENRICHED`, there were no more foreign reads there than the rest of the
  genome would predict. Most candidates fail here.
- `MIXED_TAXA`, the support comes from several organisms.
- `LOW_SUPPORT`, too few pairs to test.
- `LOW_MAPQ` and `LOW_COMPLEXITY`, the host side is not trustworthy.

With several samples, `hgt-site cohort` merges the tables and flags the
positions that turn up in nearly all of them, which is the protocol check from
above.

## What a passing site is

A lead, not a finding. It means the position has more foreign support than
chance explains, from one organism, and ideally with reads crossing the
junction. Confirm it before you believe it: PCR across the junction, or look at
the reads by eye in a genome browser.

hgt-site has been tested on simulated reads built from real genomes, including
samples where most of the chance chimeras carry the same organism as the
genuine insertion. It has not been tested on real sequencing data with a
confirmed integration, which is the one thing it still needs. If you run it on
a sample where the answer is known, a
[report](https://github.com/iraajgangavaram/hgt-site/issues/new?template=real_data_report.yml)
is genuinely useful, including one where it found nothing.

Further reading: [reading the results](interpreting.md) goes column by column,
[how it works](method.md) sets out the test, and [benchmarks](benchmarks.md)
has the measurements behind the numbers quoted here.
