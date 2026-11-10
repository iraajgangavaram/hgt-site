# What to do with the output

hgt-site gives you a table of positions. This page is about deciding which of
them are worth your time, and what to do next with the ones that are.

The short version: a site that passed is a lead, not a finding. The point of
the tool is that it has already thrown away the thousands of positions that
chance alone explains, so what is left is small enough to look at properly.

## Base composition

Three columns describe how foreign the reads at a site actually look. `gc_foreign` is their GC fraction, `gc_difference` is that minus the host's own, and `kmer_ratio` compares their 3-mer usage with the host's, divided by the difference sampling noise alone would produce at that read count.

A `kmer_ratio` near 1 means the reads are indistinguishable from host sequence. That is the useful case: it usually means the classifier called host reads foreign, which happens in repeats and with small databases, and a site built on those is suspect however good its read support looks.

**It does not separate real insertions from chance chimeras, and should not be used as a filter.** A chimera formed during library preparation carries genuine foreign sequence, so it looks exactly as foreign as an insertion does. On a benchmark run with 20 real insertions and 6 false sites, the real ones gave ratios of 2.86 to 6.72 and the false ones 2.85 to 3.94. The best single threshold that removes all the false sites also removes 19 of the 20 real ones.

Both figures are noisy at low read counts. A site with three or four reads can show a GC difference of 0.1 from chance alone, so read them alongside `n_pairs`.

The host background comes from the host-mapped mates of the anchored pairs, which are host sequence by definition. Those reads sit beside junctions rather than being spread evenly over the genome, which is a mild bias for GC and a negligible one for 3-mer usage.

## Start with the sites that passed

```bash
awk -F'\t' 'NR==1 || $NF=="PASS"' sample.hgt-site.tsv | column -ts$'\t' | less -S
```

If there are none, that is a result rather than a failure. Most host samples
contain no detectable integration, and a tool that always finds something is
not telling you anything. Check the INFO lines hgt-site printed first: if it
says few anchored pairs had a classification, the problem is the input rather
than the sample.

If there are hundreds, something is wrong with the inputs. Look at the
limitations on the [benchmarks page](benchmarks.md) and at whether the
alignment still contains the pairs where both mates mapped, which hgt-site
needs to measure fragment size.

## Reading one site

Four columns decide whether a site is interesting.

**`n_pairs` against `expected_pairs`.** This is the whole argument. The first
is how many read pairs point at this position from this taxon; the second is
how many chance alone would put there, given how common that taxon is across
the sample. A site with 20 pairs where 0.5 were expected is interesting. A
site with 4 pairs where 3 were expected is not, and will not have passed.

**`qvalue`.** The enrichment test, corrected across every position tested.
Lower is stronger. It is already filtered on, so for a site that passed, treat
it as a way to rank rather than a way to decide.

**`top_taxid_fraction`.** The share of supporting pairs from the main taxon. A
real integration gives reads from one organism. A value near 1.0 is what you
want; a lower value means several taxa point at the same position, which more
often means a repetitive or mismapping-prone region than a real event.

**`left_split_reads` and `right_split_reads`.** Reads that cross the junction
itself. Any at all is the strongest evidence in the table, because a soft
clipped read is direct sequence evidence of a join rather than an inference
from read pairs. Sites with split reads on both sides are the best candidates
you will get, and their `left_breakpoint` and `right_breakpoint` are exact to
the base.

## Sites to be sceptical of

Order the ones that passed, then look hard at these before believing them:

- **No split reads at all.** The junction is an estimate from read pairs, and
  could be a few hundred bases out. Fine for finding the site, not enough to
  design a primer around.
- **`ONE_SIDED` set.** Anchors from only one side. This is expected for a long
  insertion where the other junction is out of insert range, but it is also
  what mismapping looks like.
- **A position in a repeat or a segmental duplication.** hgt-site flags low
  complexity in the unmapped mates, not in the host region. Check the position
  in a genome browser.
- **A taxon that makes no sense for your sample.** A soil bacterium in a human
  blood sample is a reagent contaminant, and a contaminant that happens to be
  abundant can survive the test if it lands repeatedly in the same place.
- **The same position in every sample you run.** That is a property of your
  protocol or your reference, not of your biology.

## Confirming a site

Nothing in the output is confirmation. For a site you care about:

1. Look at the region in IGV or another browser, with the BED file loaded.
   Real integrations usually show a clean pile of clipped reads at one base.
2. If there are split reads, take the clipped sequence and BLAST it. It should
   match the organism the taxid says.
3. Design primers spanning the junction, using `left_breakpoint` and
   `right_breakpoint`, and run a PCR. This is the only step that settles it.

Running `--annotation` first gives you the gene at each site, which usually
decides which ones are worth that effort.

## Comparing several samples

Call each sample separately, then merge the tables:

```bash
hgt-site cohort sample*.hgt-site.tsv -o cohort.tsv
```

A position that appears in every sample, including controls, is almost
certainly an artefact of the protocol or a mismapping-prone region, however
good its statistics look in any one sample. Those come back flagged `COMMON`.
A position in one sample of a group, with split reads, is the kind of thing
worth following up. `TAXA_DIFFER` means two samples named different organisms
at the same position, which is a warning sign rather than a finding.

## When to distrust the whole run

- hgt-site reported that few anchored pairs had a classification. The read
  names in the classification file probably do not match the alignment.
- The fragment size could not be measured. The alignment has been filtered
  down to pairs with an unmapped mate, so the clustering limits are defaults
  and the junctions are rougher.
- Thousands of sites passed. Either the sample is heavily contaminated with
  one organism, or the host taxid is wrong and host reads are being counted as
  foreign.
