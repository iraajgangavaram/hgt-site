"""Scale the expected chimera count by local host coverage.

The enrichment test asks how many foreign pairs chance would put at a
position. By default it assumes they fall evenly across the genome, so the
expected count depends only on the width of the window. That is wrong
wherever host coverage is uneven, and host coverage is uneven whenever it
matters: a copy-number gain in a tumour, an amplification bias in a
low-input library, a highly expressed gene in RNA-seq.

A chimera is a host fragment joined to a foreign fragment, so the chance its
host half comes from a given region is proportional to how much host DNA from
that region is in the library. Three copies of a region contribute three times
the chimeras, and the uniform model will call the excess an insertion.

This module turns binned host coverage into a ratio of local rate to genome
mean, which multiplies the expected count. A site in a region at three times
the average coverage is then expected to hold three times the chimeras, which
is what it does hold.

The expectation is only ever raised, never lowered, so this option can cost
sensitivity but cannot cost precision. See apply_ratios for why.

Three things leave a site at a ratio of 1, tested exactly as it would be
without this option:

- An alignment filtered down to the pairs with an unmapped mate. Most of the
  host coverage has been thrown away, so what is left says nothing about it.
- Too few reads in the bins around a site to divide by.
- Local coverage at or below the genome mean.
"""

import logging

# Fully mapped pairs needed before the alignment is taken to carry host
# coverage at all.
MIN_FULLY_MAPPED = 1_000
# Mapped reads needed in the bins around a site before its own rate is used.
MIN_LOCAL_READS = 200
# A ratio this large is more likely a mapping artefact or a reference gap
# than real copy number, and it would hide a real site behind an expected
# count nothing could exceed.
MAX_RATIO = 10.0

log = logging.getLogger("hgt_site")


def unusable_reason(coverage, genome_length):
    """Why local coverage cannot be used here, or None if it can."""
    if coverage is None:
        return "the alignment was read without coverage counts"
    if coverage.mapped_reads == 0 or genome_length == 0:
        return "the alignment holds no mapped reads"
    if coverage.fully_mapped_reads < MIN_FULLY_MAPPED:
        return (
            f"only {coverage.fully_mapped_reads} fully mapped pairs are present, so the alignment "
            "appears to have been filtered down to the pairs with an unmapped mate and no longer "
            "carries the host coverage"
        )
    return None


def apply_ratios(sites, coverage, genome_length, window, contig_lengths=None):
    """Set site.depth_ratio from local host coverage. Returns how many were set.

    A site keeps a ratio of 1 where the local estimate cannot be trusted, so
    it is tested exactly as it would be without this option.
    """
    mean = coverage.mean_rate(genome_length)
    if mean <= 0:
        return 0
    scaled = 0
    for site in sites:
        # The same window the enrichment test uses, so the rate is measured
        # over the region the expected count is for.
        rate, reads = coverage.rate(
            site.chrom, site.start, site.end + window,
            (contig_lengths or {}).get(site.chrom),
        )
        if reads < MIN_LOCAL_READS:
            continue
        ratio = rate / mean
        # The bar is only ever raised, never lowered. A 10 kb bin at ordinary
        # coverage holds a few hundred reads, so the ratio carries a few per
        # cent of noise, and a ratio landing a little under 1 would lower the
        # expected count and let a marginal candidate through. On a benchmark
        # with even coverage that cost one false site in three runs. Scaling
        # up only keeps the gain where coverage is high and gives up nothing
        # where it is not: a genuinely quiet region is then tested against
        # the genome-wide rate, which is the stricter of the two.
        if ratio <= 1.0 or ratio > MAX_RATIO:
            continue
        site.depth_ratio = ratio
        scaled += 1
    return scaled
