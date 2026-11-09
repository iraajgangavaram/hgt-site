"""Regions to leave out of the results.

A protocol, a reference or a cell line will often produce the same false site
in every sample. Once a position is known to be an artefact there is no point
rediscovering it, so a BED file of such positions can be given and any
candidate overlapping one is flagged EXCLUDED rather than reported.

The reads in an excluded region still count towards the genome-wide rate used
by the enrichment test. Dropping them would lower that rate and so make every
other site look more surprising than it is, which is the opposite of what a
blacklist is for.
"""

import gzip
from bisect import bisect_left

SKIP_PREFIXES = ("#", "track", "browser")


def load_regions(path):
    """Read a BED file into {chrom: (starts, ends)} with overlaps merged.

    Only the first three columns are read. Coordinates are 0-based and
    half-open, as BED defines them and as hgt-site reports sites, so the two
    compare directly.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    by_chrom = {}
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith(SKIP_PREFIXES):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path} line {number}: a BED line needs chrom, start and end.")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(
                    f"{path} line {number}: start and end must be whole numbers."
                ) from None
            if start < 0 or end < start:
                raise ValueError(f"{path} line {number}: {start}-{end} is not a valid interval.")
            by_chrom.setdefault(fields[0], []).append((start, end))
    return {chrom: _merge(intervals) for chrom, intervals in by_chrom.items()}


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(zip(*merged))


def overlaps(regions, chrom, start, end):
    """Whether [start, end) meets any region on chrom."""
    found = regions.get(chrom)
    if not found:
        return False
    starts, ends = found
    # Merged regions are sorted and do not overlap, so their ends rise with
    # their starts. The last region beginning before `end` is therefore the
    # one reaching furthest right, and no earlier region can reach further.
    index = bisect_left(starts, end) - 1
    return index >= 0 and ends[index] > start


def flag_excluded(sites, regions):
    """Flag every site that overlaps an excluded region. Returns the count."""
    excluded = 0
    for site in sites:
        if overlaps(regions, site.chrom, site.start, site.end):
            site.flags.append("EXCLUDED")
            excluded += 1
    return excluded
