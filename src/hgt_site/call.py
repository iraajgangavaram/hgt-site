"""Run the full calling step: alignments + Kraken2 calls -> candidate site table."""

import csv
import logging

from .alignments import collect_pairs, count_by_taxon
from .kraken import TaxonFilter, load_classifications
from .sites import add_info_flags, cluster_anchors, score_sites, site_filter

log = logging.getLogger("hgt_site")

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


# Clustering limits that scale with the fragment size of the library.
SCALED_LIMITS = ("max_gap", "max_span", "window")
MIN_FRAGMENT_SAMPLE = 1000


def _percentile(values, q):
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, int(round(q / 100 * (len(ordered) - 1)))))
    return ordered[index]


def fit_limits(pairs, settings, user_set):
    """Widen the clustering limits for libraries with long fragments.

    Reads from one insertion sit within about one fragment length (less one
    read length) on each side of the junction. The defaults suit fragments of
    up to about 500 bp. If the 99th percentile of the sampled fragment sizes
    implies a wider spread, the three limits are scaled up together. They are
    never made narrower than the defaults, and limits set by the user are kept.
    """
    sizes = pairs.insert_sizes
    if len(sizes) < MIN_FRAGMENT_SAMPLE:
        log.info("Too few properly paired reads to estimate the fragment size; using the default clustering limits.")
        return settings
    fragment_99 = _percentile(sizes, 99)
    read_length = _percentile(pairs.read_lengths, 50)
    log.info(
        "Fragment size from %d pairs: median %d bp, 99th percentile %d bp; read length %d bp.",
        len(sizes), _percentile(sizes, 50), fragment_99, read_length,
    )
    needed = 2 * (fragment_99 - read_length)
    if needed <= settings["max_span"] or any(k in user_set for k in SCALED_LIMITS):
        return settings
    factor = needed / settings["max_span"]
    fitted = dict(settings)
    for key in SCALED_LIMITS:
        fitted[key] = int(round(settings[key] * factor))
    log.info(
        "Long fragments: clustering limits widened to max-gap %d, max-span %d, window %d.",
        fitted["max_gap"], fitted["max_span"], fitted["window"],
    )
    return fitted


def call_sites(alignments, kraken, host_taxid, nodes_dmp=None, exclude=(), reference=None, threads=1, **opts):
    user_set = {k for k, v in opts.items() if v is not None}
    settings = {**DEFAULTS, **{k: v for k, v in opts.items() if v is not None}}
    taxon_filter = TaxonFilter(host_taxid, nodes_dmp, exclude)
    calls = load_classifications(kraken)
    if not calls:
        log.warning("The Kraken2 file has no classified reads.")
    pairs = collect_pairs(alignments, settings["min_mapq"], reference, threads)

    with_call = [a for a in pairs.anchors if a.name in calls]
    foreign = [a for a in with_call if taxon_filter.is_foreign(calls[a.name])]
    log.info(
        "%d anchored pairs (one mate on the host, mapping quality %d or more) and %d pairs with both mates unmapped.",
        len(pairs.anchors), settings["min_mapq"], len(pairs.both_unmapped),
    )
    log.info(
        "%d anchored pairs have a Kraken2 call; %d of these are from taxa outside the host.",
        len(with_call), len(foreign),
    )
    if pairs.anchors and calls and not with_call:
        log.warning(
            "None of the anchored reads has a Kraken2 call. Check that the Kraken2 file was made "
            "from the unmapped reads of this alignment, and that the read names match."
        )
    settings = fit_limits(pairs, settings, user_set)
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


def write_bed(sites, path):
    """Write sites that passed every filter as BED, for genome browsers.

    The name column is the main taxid and the score is the number of supporting
    pairs, capped at 1000 as the BED format requires.
    """
    passed = sorted((s for s in sites if site_filter(s) == "PASS"), key=lambda s: (s.chrom, s.start))
    with open(path, "w") as handle:
        for site in passed:
            score = min(int(site.stats["n_pairs"]), 1000)
            handle.write(f"{site.chrom}\t{site.start}\t{site.end}\ttaxid_{site.stats['top_taxid']}\t{score}\t.\n")
    return len(passed)


VCF_FILTERS = {
    "LOW_SUPPORT": "Too few supporting pairs",
    "MIXED_TAXA": "Pairs from other taxa exceed what background chimeras explain",
    "LOW_MAPQ": "Low mean mapping quality of the anchors",
    "LOW_COMPLEXITY": "Low-complexity sequence in the unmapped mates",
    "NOT_ENRICHED": "Not enriched over chance chimeras (q-value above the threshold)",
}


def _junction(site):
    """1-based position of the estimated junction, falling back to the site start."""
    left, right = site.stats.get("left_junction"), site.stats.get("right_junction")
    if left != "" and left is not None:
        return int(left)
    if right != "" and right is not None:
        return int(right) + 1
    return site.start + 1


def write_vcf(sites, path, version="unknown"):
    """Write every candidate site as a VCF insertion record.

    FILTER is PASS for sites that passed every filter, otherwise the failing
    flags. The reference base is not known to hgt-site, so REF is N and ALT is
    the symbolic allele <INS>.
    """
    ordered = sorted(sites, key=lambda s: (s.chrom, _junction(s)))
    contigs = sorted({s.chrom for s in ordered})
    with open(path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##source=hgt-site {version}\n")
        for contig in contigs:
            handle.write(f"##contig=<ID={contig}>\n")
        handle.write('##ALT=<ID=INS,Description="Insertion of foreign sequence">\n')
        for name, text in VCF_FILTERS.items():
            handle.write(f'##FILTER=<ID={name},Description="{text}">\n')
        for line in (
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
            '##INFO=<ID=SITE_START,Number=1,Type=Integer,Description="Start of the host region covered by the anchors (1-based)">',
            '##INFO=<ID=SITE_END,Number=1,Type=Integer,Description="End of the host region covered by the anchors">',
            '##INFO=<ID=TAXID,Number=1,Type=Integer,Description="NCBI taxid of the main foreign taxon">',
            '##INFO=<ID=NPAIRS,Number=1,Type=Integer,Description="Number of supporting read pairs">',
            '##INFO=<ID=NFWD,Number=1,Type=Integer,Description="Supporting pairs with a forward anchor">',
            '##INFO=<ID=NREV,Number=1,Type=Integer,Description="Supporting pairs with a reverse anchor">',
            '##INFO=<ID=QVAL,Number=1,Type=Float,Description="Benjamini-Hochberg q-value of the enrichment test">',
            '##INFO=<ID=ONE_SIDED,Number=0,Type=Flag,Description="Anchors from only one side of the junction">',
        ):
            handle.write(line + "\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, site in enumerate(ordered, start=1):
            st = site.stats
            pos = _junction(site)
            info = [
                "SVTYPE=INS", f"END={pos}", f"SITE_START={site.start + 1}", f"SITE_END={site.end}",
                f"TAXID={st['top_taxid']}", f"NPAIRS={st['n_pairs']}", f"NFWD={st['n_forward']}",
                f"NREV={st['n_reverse']}", f"QVAL={st['qvalue']:.4g}",
            ]
            if "ONE_SIDED" in site.flags:
                info.append("ONE_SIDED")
            handle.write(
                f"{site.chrom}\t{pos}\thgtsite_{i}\tN\t<INS>\t.\t{site_filter(site)}\t{';'.join(info)}\n"
            )
    return len(ordered)
