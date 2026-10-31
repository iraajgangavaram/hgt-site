"""Run the full calling step: alignments + Kraken2 calls -> candidate site table."""

import csv

from .alignments import collect_pairs, count_by_taxon
from .kraken import TaxonFilter, load_classifications
from .sites import add_info_flags, cluster_anchors, score_sites, site_filter

DEFAULTS = {
    "min_mapq": 20,
    "max_gap": 300,
    "max_span": 800,
    "window": 600,
    "min_support": 3,
    "min_taxon_fraction": 0.8,
    "max_other_pvalue": 0.01,
    "min_mean_mapq": 20,
    "min_complexity": 0.5,
    "max_qvalue": 0.05,
}

COLUMNS = [
    "chrom",
    "start",
    "end",
    "left_junction",
    "right_junction",
    "n_pairs",
    "n_forward",
    "n_reverse",
    "top_taxid",
    "top_taxid_fraction",
    "n_taxa",
    "mean_mapq",
    "mean_mate_complexity",
    "taxon_anchored_pairs",
    "taxon_unanchored_pairs",
    "expected_pairs",
    "pvalue",
    "qvalue",
    "expected_other_pairs",
    "other_pvalue",
    "flags",
    "filter",
]


def call_sites(alignments, kraken, host_taxid, nodes_dmp=None, exclude=(), reference=None, **opts):
    settings = {**DEFAULTS, **{k: v for k, v in opts.items() if v is not None}}
    taxon_filter = TaxonFilter(host_taxid, nodes_dmp, exclude)
    calls = load_classifications(kraken)
    pairs = collect_pairs(alignments, settings["min_mapq"], reference)

    foreign = [
        a
        for a in pairs.anchors
        if a.name in calls and taxon_filter.is_foreign(calls[a.name])
    ]
    if not foreign:
        return []

    background = count_by_taxon(pairs.both_unmapped, calls, taxon_filter)
    genome_length = sum(pairs.contig_lengths.values())
    if genome_length == 0:
        raise ValueError("The alignment header lists no reference contigs and lengths.")

    sites = cluster_anchors(foreign, settings["max_gap"], settings["max_span"])
    sites = score_sites(
        sites, calls, pairs.unmapped_mate_seqs, background, genome_length,
        settings["window"], settings,
    )
    return add_info_flags(sites)


def _fmt(value):
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value)


def write_table(sites, path):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for site in sorted(sites, key=lambda s: (s.chrom, s.start)):
            row = {"chrom": site.chrom, "start": site.start, "end": site.end, **site.stats}
            row["flags"] = ";".join(site.flags) if site.flags else "."
            row["filter"] = site_filter(site)
            writer.writerow([_fmt(row[c]) for c in COLUMNS])
