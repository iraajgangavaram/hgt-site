"""Group anchored reads into candidate sites, score them and flag likely artefacts."""

import math
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import product

KMER = 3
_ALL_KMERS = len(list(product("ACGT", repeat=KMER)))


def sequence_complexity(seq, k=KMER):
    """Shannon entropy of the k-mer counts, scaled to 0-1.

    Low values point to repeats or low-complexity sequence, which Kraken2
    tends to misclassify.
    """
    seq = seq.upper()
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    kmers = [km for km in kmers if set(km) <= set("ACGT")]
    if len(kmers) < 2:
        return 0.0
    counts = Counter(kmers)
    total = len(kmers)
    entropy = -sum((c / total) * math.log(c / total) for c in counts.values())
    max_entropy = math.log(min(_ALL_KMERS, total))
    return entropy / max_entropy if max_entropy > 0 else 0.0


def poisson_sf(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    below = 0.0
    for i in range(k):
        below += math.exp(i * math.log(lam) - lam - math.lgamma(i + 1))
    return max(0.0, 1.0 - below)


def benjamini_hochberg(pvalues, n_tests=None):
    """Benjamini-Hochberg q-values.

    n_tests can be set higher than the number of p-values when more tests
    were implicitly run than were reported.
    """
    n = max(len(pvalues), n_tests or 0)
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    qvalues = [0.0] * len(pvalues)
    running = 1.0
    for rank in range(len(pvalues), 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * n / rank)
        qvalues[i] = running
    return qvalues


@dataclass
class Site:
    chrom: str
    start: int
    end: int
    anchors: list = field(default_factory=list)
    taxa: Counter = field(default_factory=Counter)
    complexity: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)
    flags: list = field(default_factory=list)


def cluster_anchors(anchors, max_gap, max_span):
    """Cluster anchors by the side that faces the junction.

    Anchors join a cluster while they are within max_gap of the previous
    anchor and within max_span of the first one. The span limit stops
    scattered chimeric reads chaining into one long false cluster, since
    real anchors sit within about one fragment length of the junction.
    """
    by_chrom = defaultdict(list)
    for anchor in anchors:
        by_chrom[anchor.chrom].append(anchor)
    sites = []
    for chrom in sorted(by_chrom):
        # Sorting on the whole anchor, not only the junction side, so that
        # anchors that face the junction at the same base always fall in the
        # same order. Otherwise clusters could depend on the order of the
        # records in the file, and a sorted BAM would give slightly different
        # sites from the same reads in a SAM file.
        group = sorted(
            by_chrom[chrom], key=lambda a: (a.junction_side, a.start, a.end, a.reverse, a.name)
        )
        current = [group[0]]
        for anchor in group[1:]:
            close = anchor.junction_side - current[-1].junction_side <= max_gap
            within = anchor.junction_side - current[0].junction_side <= max_span
            if close and within:
                current.append(anchor)
            else:
                sites.append(current)
                current = [anchor]
        sites.append(current)
    return [
        Site(chrom=c[0].chrom, start=min(a.start for a in c), end=max(a.end for a in c), anchors=c)
        for c in sites
    ]


def score_sites(sites, calls, mate_seqs, background, genome_length, window, thresholds):
    """Fill in the statistics and flags for each site.

    The enrichment test asks whether a site has more anchored pairs for its
    main taxon than we would expect if that taxon's anchored pairs were
    spread at random across the host genome, as chimeras made during library
    preparation would be.
    """
    anchored_per_taxon = Counter()
    for site in sites:
        for anchor in site.anchors:
            anchored_per_taxon[calls[anchor.name]] += 1

    pvalues = []
    for site in sites:
        for anchor in site.anchors:
            site.taxa[calls[anchor.name]] += 1
            seq = mate_seqs.get(anchor.name)
            if seq:
                site.complexity.append(sequence_complexity(seq))
        top_taxid, top_count = site.taxa.most_common(1)[0]
        n_pairs = len(site.anchors)
        n_rev = sum(a.reverse for a in site.anchors)
        n_fwd = n_pairs - n_rev
        fwd_ends = [a.end for a in site.anchors if not a.reverse]
        rev_starts = [a.start for a in site.anchors if a.reverse]
        span = (site.end - site.start) + window
        expected = anchored_per_taxon[top_taxid] * span / genome_length
        pvalue = poisson_sf(top_count, expected)
        pvalues.append(pvalue)
        # Pairs from other taxa that chance chimeras alone would put in this window.
        other_total = sum(anchored_per_taxon.values()) - anchored_per_taxon[top_taxid]
        expected_other = other_total * span / genome_length
        other_count = n_pairs - top_count
        site.stats = {
            "n_pairs": n_pairs,
            "n_forward": n_fwd,
            "n_reverse": n_rev,
            "left_junction": max(fwd_ends) if fwd_ends else "",
            "right_junction": min(rev_starts) if rev_starts else "",
            "top_taxid": top_taxid,
            "top_taxid_fraction": top_count / n_pairs,
            "n_taxa": len(site.taxa),
            "mean_mapq": sum(a.mapq for a in site.anchors) / n_pairs,
            "mean_mate_complexity": (
                sum(site.complexity) / len(site.complexity) if site.complexity else ""
            ),
            "taxon_anchored_pairs": anchored_per_taxon[top_taxid],
            "taxon_unanchored_pairs": background.get(top_taxid, 0),
            "expected_pairs": expected,
            "pvalue": pvalue,
            "expected_other_pairs": expected_other,
            "other_pvalue": poisson_sf(other_count, expected_other),
        }

    # Every window of the genome could have produced a cluster for every
    # taxon with anchored reads, not only the windows where one was seen.
    # Counting only observed clusters would make the correction too lenient.
    n_windows = math.ceil(genome_length / window)
    n_tests = n_windows * len(anchored_per_taxon)
    for site, qvalue in zip(sites, benjamini_hochberg(pvalues, n_tests)):
        site.stats["qvalue"] = qvalue
        s, t = site.stats, thresholds
        if s["n_pairs"] < t["min_support"]:
            site.flags.append("LOW_SUPPORT")
        # A site is only mixed if it holds more reads from other taxa than
        # background chimeras explain. Without this, real sites in noisy
        # samples fail because a few random chimeras land inside them.
        excess_other = s["other_pvalue"] < t["max_other_pvalue"]
        if s["top_taxid_fraction"] < t["min_taxon_fraction"] and excess_other:
            site.flags.append("MIXED_TAXA")
        if s["mean_mapq"] < t["min_mean_mapq"]:
            site.flags.append("LOW_MAPQ")
        if s["mean_mate_complexity"] != "" and s["mean_mate_complexity"] < t["min_complexity"]:
            site.flags.append("LOW_COMPLEXITY")
        if qvalue > t["max_qvalue"]:
            site.flags.append("NOT_ENRICHED")
    return sites


# Flags that make a site fail. ONE_SIDED is reported for information only,
# because a long insertion can leave only one junction within insert range.
INFO_FLAGS = {"ONE_SIDED"}


def add_info_flags(sites):
    for site in sites:
        if site.stats["n_forward"] == 0 or site.stats["n_reverse"] == 0:
            site.flags.append("ONE_SIDED")
    return sites


def site_filter(site):
    failing = [f for f in site.flags if f not in INFO_FLAGS]
    return "PASS" if not failing else ";".join(failing)
