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

Two things make the estimate untrustworthy, and in both the ratio falls back
to 1, which is the uniform model:

- An alignment filtered down to the pairs with an unmapped mate. Most of the
  host coverage has been thrown away, so what is left says nothing about it.
- Too few reads in the bins around a site to divide by.
"""

import logging

# Fully mapped pairs needed before the alignment is taken to carry host
# coverage at all.
MIN_FULLY_MAPPED = 1_000
# Mapped reads needed in the bins around a site before its own rate is used.
MIN_LOCAL_READS = 200
# Ratios outside this range are more likely a mapping artefact, a reference
# gap or a contig edge than real copy number, and a huge ratio would hide a
# real site behind an expected count nothing could exceed.
MIN_RATIO, MAX_RATIO = 0.1, 10.0

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
        if not MIN_RATIO <= ratio <= MAX_RATIO:
            continue
        site.depth_ratio = ratio
        scaled += 1
    return scaled
