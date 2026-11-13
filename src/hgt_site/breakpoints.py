"""Find the exact junction of each candidate site from split reads.

Read pairs place a junction only to within about one fragment length. Reads
that cross the junction itself align to the host up to the junction and are
soft-clipped beyond it, so many clips ending at the same base give the
junction to the base pair.

On the side covered by forward anchors the host sequence ends at the
junction, so crossing reads are clipped on their right. On the side covered
by reverse anchors the host sequence starts at the junction, so crossing
reads are clipped on their left. The two positions can differ by a few bases
when the insertion duplicated or deleted host sequence.

The BAM or CRAM file must be indexed, so that only the reads around each
site are read. Plain SAM files are scanned in full, which is only sensible
for small files.
"""

import logging
from bisect import bisect_left
from collections import Counter

log = logging.getLogger("hgt_site")

# Shortest soft clip that counts as crossing a junction. Shorter clips are
# mostly sequencing errors and adapter remains near read ends.
MIN_CLIP = 20
# Minimum number of reads that must clip at the same base.
MIN_SPLIT_READS = 2
# How far a clip may sit inside the anchors, to allow for anchors that are
# themselves clipped at the junction.
OVERLAP = 10

SOFT_CLIP, HARD_CLIP = 4, 5
UNMAPPED, SECONDARY, QCFAIL, DUPLICATE = 0x4, 0x100, 0x200, 0x400
SKIP = UNMAPPED | SECONDARY | QCFAIL | DUPLICATE


def _search_windows(site, window):
    """Host ranges where crossing reads could end, one per side of the junction."""
    left = site.stats.get("left_junction")
    right = site.stats.get("right_junction")
    windows = {}
    if left not in ("", None):
        windows["left"] = (int(left) - OVERLAP, int(left) + window)
    if right not in ("", None):
        windows["right"] = (max(0, int(right) - window), int(right) + OVERLAP)
    return windows


def _clips(cigar):
    """Lengths of the clips at the start and end of an alignment."""
    if not cigar:
        return 0, 0
    first, last = cigar[0], cigar[-1]
    start = first[1] if first[0] in (SOFT_CLIP, HARD_CLIP) else 0
    end = last[1] if last[0] in (SOFT_CLIP, HARD_CLIP) and len(cigar) > 1 else 0
    return start, end


def _count(reads, windows, min_mapq):
    """Count right clips in the left window and left clips in the right window."""
    counts = {"left": Counter(), "right": Counter()}
    for _chrom, flag, mapq, start, end, cigar in reads:
        if flag & SKIP or mapq < min_mapq:
            continue
        clip_start, clip_end = _clips(cigar)
        if clip_end >= MIN_CLIP and "left" in windows:
            lo, hi = windows["left"]
            if lo <= end <= hi:
                counts["left"][end] += 1
        if clip_start >= MIN_CLIP and "right" in windows:
            lo, hi = windows["right"]
            if lo <= start <= hi:
                counts["right"][start] += 1
    return counts


def _best(counter, anchor_edge):
    """Most common clip position, ties going to the one nearest the anchors."""
    if not counter:
        return "", 0
    position, n = max(counter.items(), key=lambda kv: (kv[1], -abs(kv[0] - anchor_edge)))
    if n < MIN_SPLIT_READS:
        return "", 0
    return position, n


def _record(site, counts):
    left_edge = site.stats.get("left_junction")
    right_edge = site.stats.get("right_junction")
    pos, n = _best(counts["left"], int(left_edge) if left_edge not in ("", None) else 0)
    site.stats["left_breakpoint"], site.stats["left_split_reads"] = pos, n
    pos, n = _best(counts["right"], int(right_edge) if right_edge not in ("", None) else 0)
    site.stats["right_breakpoint"], site.stats["right_split_reads"] = pos, n


def _parse_cigar(text):
    ops = {"M": 0, "I": 1, "D": 2, "N": 3, "S": 4, "H": 5, "P": 6, "=": 7, "X": 8}
    cigar, number = [], ""
    for char in text:
        if char.isdigit():
            number += char
        else:
            cigar.append((ops[char], int(number)))
            number = ""
    return cigar


def _sam_clipped_reads(path, min_mapq):
    """Records of a SAM file that are clipped enough to cross a junction.

    The cheap text checks come first: a clipped alignment must have an S or H
    in its CIGAR, and most reads have neither, so very few lines are parsed.
    """
    with open(path) as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.split("\t", 6)
            cigar_text = cols[5]
            if "S" not in cigar_text and "H" not in cigar_text:
                continue
            flag = int(cols[1])
            if flag & SKIP or int(cols[4]) < min_mapq:
                continue
            cigar = _parse_cigar(cigar_text)
            start = int(cols[3]) - 1
            end = start + sum(n for op, n in cigar if op in (0, 2, 3, 7, 8))
            yield cols[2], flag, int(cols[4]), start, end, cigar


def _empty(site):
    site.stats.update(left_breakpoint="", left_split_reads=0, right_breakpoint="", right_split_reads=0)


def find_breakpoints(
    sites, path, window, min_mapq=20, reference=None, threads=1, targets=None, index=None
):
    """Add left/right breakpoints and split-read counts to each site's stats.

    Every site gets the columns, but only the sites in targets are searched.
    A noisy sample can hold thousands of clusters that fail their other
    filters, and searching those would cost time for results nobody uses.
    Returns True if the alignments could be searched.
    """
    for site in sites:
        _empty(site)
    targets = list(sites if targets is None else targets)
    if not targets:
        return True
    path = str(path)

    if path.endswith(".sam"):
        # A SAM file cannot be indexed, so it is scanned once. Sites are held
        # in a position index, so each clipped read is compared with the few
        # sites near it rather than with all of them.
        index = {}
        windows_by_site = {}
        for site in targets:
            windows = _search_windows(site, window)
            windows_by_site[id(site)] = windows
            for lo, hi in windows.values():
                index.setdefault(site.chrom, []).append((lo, hi, site))
        for rows in index.values():
            rows.sort(key=lambda row: row[0])
        starts = {chrom: [row[0] for row in rows] for chrom, rows in index.items()}
        widest = max((hi - lo for rows in index.values() for lo, hi, _ in rows), default=0)

        counts = {id(s): {"left": Counter(), "right": Counter()} for s in targets}
        for read in _sam_clipped_reads(path, min_mapq):
            rows = index.get(read[0])
            if not rows:
                continue
            chrom_starts = starts[read[0]]
            first = bisect_left(chrom_starts, read[3] - widest)
            seen = set()
            for lo, hi, site in rows[first:]:
                if lo > read[4]:
                    break
                if hi < read[3] or id(site) in seen:
                    continue
                seen.add(id(site))
                found = _count([read], windows_by_site[id(site)], min_mapq)
                for side in ("left", "right"):
                    counts[id(site)][side].update(found[side])
        for site in targets:
            _record(site, counts[id(site)])
        return True

    import pysam

    mode = "rc" if path.endswith(".cram") else "rb"
    with pysam.AlignmentFile(
        path, mode, reference_filename=reference, threads=threads, index_filename=index
    ) as bam:
        if not bam.has_index():
            log.info(
                "The alignment file has no index, so junctions are placed from read pairs only. "
                "Index it (samtools index) to find exact junctions from split reads."
            )
            return False
        for site in targets:
            windows = _search_windows(site, window)
            if not windows:
                continue
            lo = min(w[0] for w in windows.values())
            hi = max(w[1] for w in windows.values())
            reads = (
                (r.reference_name, r.flag, r.mapping_quality, r.reference_start, r.reference_end, r.cigartuples)
                for r in bam.fetch(site.chrom, max(0, lo - 1), hi + 1)
            )
            _record(site, _count(reads, windows, min_mapq))
    return True
