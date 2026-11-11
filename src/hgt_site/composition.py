"""Base composition of the foreign reads at a site, against the host's own.

Sequence from a different organism usually differs from the host in GC content
and in short k-mer usage, which is one of the oldest signals for horizontal
transfer. It is reported here rather than used as a filter, for two reasons.

A real donor can happen to share the host's composition, so a small difference
is not evidence against a site. And a large difference is not evidence for
one either, because a chance chimera carries genuine foreign sequence too and
so looks exactly as foreign as a real insertion does.

What the numbers are good for is the opposite case: foreign reads whose
composition is indistinguishable from the host are often host reads the
classifier got wrong, particularly in repeats. A site built on those is
suspect however good its read support looks.

The host background is taken from the host-mapped mates of the anchored pairs,
which are host sequence by definition and are already being read, so nothing
extra is scanned. Those reads sit beside junctions rather than being spread
evenly over the genome, which is a mild bias for GC and a negligible one for
3-mer usage.
"""

import math
from collections import Counter

K = 3
# Enough reads to settle a GC fraction and a 64-way 3-mer profile. Thousands
# of reads put the standard error on GC well under a percentage point.
SAMPLE_CAP = 20_000


class Background:
    """GC and k-mer usage accumulated over a sample of host reads."""

    __slots__ = ("gc", "bases", "kmers", "reads", "cap")

    def __init__(self, cap=SAMPLE_CAP):
        self.gc = 0
        self.bases = 0
        self.kmers = Counter()
        self.reads = 0
        self.cap = cap

    def add(self, seq):
        if self.reads >= self.cap or not seq:
            return
        self.reads += 1
        seq = seq.upper()
        self.gc += seq.count("G") + seq.count("C")
        self.bases += len(seq) - seq.count("N")
        self.kmers.update(kmers(seq))

    @property
    def gc_fraction(self):
        return self.gc / self.bases if self.bases else None

    @property
    def profile(self):
        return normalise(self.kmers)


def kmers(seq, k=K):
    """The k-mers of a sequence, skipping any containing a base other than ACGT."""
    allowed = set("ACGT")
    return [
        seq[i : i + k] for i in range(len(seq) - k + 1) if set(seq[i : i + k]) <= allowed
    ]


def normalise(counts):
    """Counts as fractions of their total."""
    total = sum(counts.values())
    return {kmer: count / total for kmer, count in counts.items()} if total else {}


def gc_fraction(seqs):
    """GC fraction over several sequences, ignoring Ns."""
    gc = bases = 0
    for seq in seqs:
        seq = seq.upper()
        gc += seq.count("G") + seq.count("C")
        bases += len(seq) - seq.count("N")
    return gc / bases if bases else None


def distance(first, second):
    """Total variation distance between two k-mer profiles, 0 to 1.

    Half the sum of absolute differences, so two identical profiles give 0 and
    two with no k-mer in common give 1. It needs no pseudocounts and does not
    blow up when a k-mer is missing from one side, which both happen with a
    few hundred reads.
    """
    if not first or not second:
        return None
    # Summed in sorted order so the result does not depend on how the sets
    # happen to iterate, which would otherwise make the last digits differ
    # between runs and between platforms.
    return sum(
        abs(first.get(k, 0.0) - second.get(k, 0.0)) for k in sorted(set(first) | set(second))
    ) / 2


def expected_noise(profile, n):
    """How far a sample of n k-mers drawn from `profile` sits from it by chance.

    A site has a few hundred to a few thousand k-mers spread over 64 bins, so
    an appreciable part of any distance is sampling noise, and the amount
    depends on how many reads the site happens to have. Dividing the measured
    distance by this makes the number comparable between sites.

    Each bin's count is close to normal, so the mean absolute deviation is
    sqrt(2 p (1 - p) / (pi n)), and the total variation distance is half the
    sum of those.
    """
    if not profile or n <= 0:
        return None
    return 0.5 * sum(
        math.sqrt(2 * p * (1 - p) / (math.pi * n)) for p in profile.values() if p > 0
    )


def summarise(seqs, host):
    """GC and k-mer composition for a site's foreign reads against the host.

    Returns (gc, gc_difference, kmer_ratio), where kmer_ratio is the k-mer
    distance over the distance expected from sampling noise alone. Around 1
    means the reads are indistinguishable from host sequence; well above 1
    means their composition genuinely differs. None where there is nothing to
    compare.
    """
    seqs = [s for s in seqs if s]
    if not seqs:
        return None, None, None
    gc = gc_fraction(seqs)
    host_gc = host.gc_fraction if host else None
    difference = None if gc is None or host_gc is None else gc - host_gc
    if not host:
        return gc, difference, None
    counts = Counter()
    for seq in seqs:
        counts.update(kmers(seq.upper()))
    total = sum(counts.values())
    host_profile = host.profile
    measured = distance(normalise(counts), host_profile)
    noise = expected_noise(host_profile, total)
    ratio = None if not measured or not noise else measured / noise
    # Rounded because these are descriptive numbers read by a person, and
    # carrying sixteen digits only invites differences between platforms.
    return (
        None if gc is None else round(gc, 4),
        None if difference is None else round(difference, 4),
        None if ratio is None else round(ratio, 3),
    )
