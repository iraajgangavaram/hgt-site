"""Run the full calling step: alignments + Kraken2 calls -> candidate site table."""

import csv
import logging
import math
import time

from .alignments import collect_pairs, count_by_taxon
from .breakpoints import find_breakpoints
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
    "left_breakpoint",
    "left_split_reads",
    "right_breakpoint",
    "right_split_reads",
    "n_pairs",
    "n_forward",
    "n_reverse",
    "gene",
    "gene_id",
    "gene_distance",
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


def call_sites(
    alignments, kraken, host_taxid, nodes_dmp=None, exclude=(), reference=None, threads=1, split_reads=True,
    annotation=None, max_gene_distance=100_000, classifier="auto", summary=None, **opts
):
    """Call candidate sites. If summary is a dict, read counts are added to it."""
    user_set = {k for k, v in opts.items() if v is not None}
    settings = {**DEFAULTS, **{k: v for k, v in opts.items() if v is not None}}
    taxon_filter = TaxonFilter(host_taxid, nodes_dmp, exclude)
    calls = load_classifications(kraken, classifier)
    if not calls:
        log.warning("The classification file has no classified reads.")
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
    if summary is not None:
        summary.update(
            anchored_pairs=len(pairs.anchors), anchored_pairs_with_call=len(with_call),
            foreign_anchored_pairs=len(foreign), both_unmapped_pairs=len(pairs.both_unmapped),
            contig_lengths=dict(pairs.contig_lengths),
        )
    settings = fit_limits(pairs, settings, user_set)
    if summary is not None:
        # After fit_limits, so these are the limits actually used rather than
        # the ones asked for. The VCF header records them.
        summary["settings"] = dict(settings)
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
    if split_reads:
        # Only the sites that passed their other filters are worth searching.
        passing = [s for s in sites if site_filter(s) == "PASS"]
        find_breakpoints(
            sites, alignments, settings["window"], settings["min_mapq"], reference, threads, targets=passing
        )
        found = sum(1 for s in passing if s.stats["left_split_reads"] or s.stats["right_split_reads"])
        log.info("Exact junctions from split reads found for %d of %d sites that passed.", found, len(passing))
    else:
        for site in sites:
            site.stats.update(left_breakpoint="", left_split_reads=0, right_breakpoint="", right_split_reads=0)
    sites = add_info_flags(sites)
    if annotation:
        from .annotate import annotate

        found = annotate(sites, annotation, _junction, max_gene_distance)
        log.info("A gene was found at or near %d of %d candidate sites.", found, len(sites))
    else:
        for site in sites:
            site.stats.setdefault("gene", "")
            site.stats.setdefault("gene_id", "")
            site.stats.setdefault("gene_distance", "")
    return sites


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


def _given(value):
    return value not in ("", None)


def _junction(site):
    """1-based position of the last host base before the insertion.

    Split-read breakpoints are used when there are any, then the edges of the
    anchors, then the site start. All of these are 0-based boundaries, which
    equal the 1-based position of the base before the boundary.
    """
    st = site.stats
    for key in ("left_breakpoint", "right_breakpoint", "left_junction", "right_junction"):
        if _given(st.get(key)):
            return max(1, int(st[key]))
    return site.start + 1


def _precise(site):
    st = site.stats
    return bool(st.get("left_split_reads") or st.get("right_split_reads"))


# The highest QUAL we report. A q-value of zero after rounding would otherwise
# give infinity, and no caller gains anything from a larger number.
MAX_QUAL = 1000


def _qual(site):
    """Phred-scale the q-value, so the VCF can be filtered on QUAL.

    A site that fails the enrichment test still gets its score, because FILTER
    already says it failed and a reader may want to rank the failures.
    """
    q = site.stats.get("qvalue")
    if q is None or q == "":
        return "."
    q = float(q)
    if q <= 0:
        return str(MAX_QUAL)
    if q >= 1:
        return "0"
    return f"{min(MAX_QUAL, -10 * math.log10(q)):.1f}"


def _confidence_interval(site, pos):
    """How far the junction could be from pos, as the VCF CIPOS pair.

    With split reads the position is exact and there is no interval. Without
    them, the two sides of the junction give two estimates, and the junction
    lies between them. If only one side has anchors there is no second
    estimate, so the interval is the whole region the anchors cover, which is
    the most that can honestly be said.
    """
    st = site.stats
    estimates = [int(st[key]) for key in ("left_junction", "right_junction") if _given(st.get(key))]
    if len(estimates) < 2:
        estimates = [site.start + 1, site.end]
    low = min(min(estimates), pos) - pos
    high = max(max(estimates), pos) - pos
    return low, high


def _settings_lines(settings):
    for key in sorted(settings):
        yield f"##hgt_siteSetting=<ID={key},Value={settings[key]}>"


def write_vcf(sites, path, version="unknown", contigs=None, command=None, settings=None,
              sample=None, reference=None):
    """Write every candidate site as a VCF insertion record.

    FILTER is PASS for sites that passed every filter, otherwise the failing
    flags. The reference base is not known to hgt-site, so REF is N and ALT is
    the symbolic allele <INS>.

    contigs maps contig name to length, taken from the alignment header, and is
    written into the header so other tools can index the file. command and
    settings are recorded so a file can be traced back to how it was made. If
    path ends in .gz the output is block compressed and indexed with tabix.

    A single sample column is written, with the paired and split read counts.
    The genotype is always ./. because hgt-site counts the reads supporting a
    junction and does not call how many copies of the host carry it.
    """
    ordered = sorted(sites, key=lambda s: (s.chrom, _junction(s)))
    contigs = contigs or {}
    names = sorted({s.chrom for s in ordered} | set(contigs))
    sample = sample or "SAMPLE"
    compress = str(path).endswith(".gz")
    text_path = str(path)[: -len(".gz")] if compress else str(path)

    with open(text_path, "w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##fileDate={time.strftime('%Y%m%d')}\n")
        handle.write(f"##source=hgt-site {version}\n")
        if reference:
            handle.write(f"##reference={reference}\n")
        if command:
            handle.write(f"##hgt_siteCommand={command}\n")
        if settings:
            for line in _settings_lines(settings):
                handle.write(line + "\n")
        for name in names:
            length = contigs.get(name)
            if length:
                handle.write(f"##contig=<ID={name},length={length}>\n")
            else:
                handle.write(f"##contig=<ID={name}>\n")
        handle.write('##ALT=<ID=INS,Description="Insertion of foreign sequence">\n')
        for name, text in VCF_FILTERS.items():
            handle.write(f'##FILTER=<ID={name},Description="{text}">\n')
        for line in (
            '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
            '##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">',
            '##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">',
            '##INFO=<ID=SITE_START,Number=1,Type=Integer,Description="Start of the host region covered by the anchors (1-based)">',
            '##INFO=<ID=SITE_END,Number=1,Type=Integer,Description="End of the host region covered by the anchors">',
            '##INFO=<ID=TAXID,Number=1,Type=Integer,Description="NCBI taxid of the main foreign taxon">',
            '##INFO=<ID=NPAIRS,Number=1,Type=Integer,Description="Number of supporting read pairs">',
            '##INFO=<ID=NFWD,Number=1,Type=Integer,Description="Supporting pairs with a forward anchor">',
            '##INFO=<ID=NREV,Number=1,Type=Integer,Description="Supporting pairs with a reverse anchor">',
            '##INFO=<ID=QVAL,Number=1,Type=Float,Description="Benjamini-Hochberg q-value of the enrichment test">',
            '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene at the insertion point, or the nearest gene within 100 kb">',
            '##INFO=<ID=GENEDIST,Number=1,Type=Integer,Description="Distance in bp from the insertion point to that gene, 0 if inside it">',
            '##INFO=<ID=SPLIT_LEFT,Number=1,Type=Integer,Description="Split reads at the junction on the forward-anchor side">',
            '##INFO=<ID=SPLIT_RIGHT,Number=1,Type=Integer,Description="Split reads at the junction on the reverse-anchor side">',
            '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="No split reads; position estimated from read pairs">',
            '##INFO=<ID=ONE_SIDED,Number=0,Type=Flag,Description="Anchors from only one side of the junction">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype, always ./. because hgt-site does not genotype insertions">',
            '##FORMAT=<ID=PR,Number=1,Type=Integer,Description="Read pairs supporting the junction">',
            '##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Split reads crossing the junction">',
        ):
            handle.write(line + "\n")
        handle.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
        for i, site in enumerate(ordered, start=1):
            st = site.stats
            pos = _junction(site)
            split = int(st.get("left_split_reads") or 0) + int(st.get("right_split_reads") or 0)
            info = [
                "SVTYPE=INS", f"END={pos}", f"SITE_START={site.start + 1}", f"SITE_END={site.end}",
                f"TAXID={st['top_taxid']}", f"NPAIRS={st['n_pairs']}", f"NFWD={st['n_forward']}",
                f"NREV={st['n_reverse']}", f"QVAL={st['qvalue']:.4g}",
                f"SPLIT_LEFT={st.get('left_split_reads', 0)}", f"SPLIT_RIGHT={st.get('right_split_reads', 0)}",
            ]
            if st.get("gene"):
                info.append(f"GENE={st['gene']}")
                info.append(f"GENEDIST={st['gene_distance']}")
            if not _precise(site):
                low, high = _confidence_interval(site, pos)
                info.append(f"CIPOS={low},{high}")
                info.append(f"CIEND={low},{high}")
                info.append("IMPRECISE")
            if "ONE_SIDED" in site.flags:
                info.append("ONE_SIDED")
            handle.write(
                f"{site.chrom}\t{pos}\thgtsite_{i}\tN\t<INS>\t{_qual(site)}\t{site_filter(site)}\t"
                f"{';'.join(info)}\tGT:PR:SR\t./.:{st['n_pairs']}:{split}\n"
            )

    if compress:
        _compress_vcf(text_path, str(path))
    return len(ordered)


def _compress_vcf(text_path, gz_path):
    """Block compress a VCF and index it, so pipelines can use it directly."""
    import os

    import pysam

    pysam.tabix_compress(text_path, gz_path, force=True)
    pysam.tabix_index(gz_path, preset="vcf", force=True)
    os.remove(text_path)


def write_multiqc(sites, path, sample, summary=None):
    """Write a MultiQC custom-content table with one row for this sample.

    MultiQC picks the file up when its name ends in _mqc.json.
    """
    import json

    passed = [s for s in sites if site_filter(s) == "PASS"]
    row = dict(summary or {})
    row.update(
        candidate_sites=len(sites),
        passed_sites=len(passed),
        passed_with_exact_junction=sum(1 for s in passed if _precise(s)),
    )
    headers = {
        "passed_sites": {"title": "HGT sites", "description": "Candidate sites that passed every filter", "format": "{:,.0f}"},
        "passed_with_exact_junction": {"title": "Exact junctions", "description": "Passed sites with split reads at the junction", "format": "{:,.0f}"},
        "candidate_sites": {"title": "Candidates", "description": "All candidate sites, before filtering", "format": "{:,.0f}"},
        "foreign_anchored_pairs": {"title": "Foreign anchored pairs", "description": "Pairs with one mate on the host and the other classified outside the host", "format": "{:,.0f}"},
        "anchored_pairs": {"title": "Anchored pairs", "description": "Pairs with one mate on the host and the other unmapped", "format": "{:,.0f}", "hidden": True},
        "anchored_pairs_with_call": {"title": "Anchored pairs with a Kraken2 call", "format": "{:,.0f}", "hidden": True},
        "both_unmapped_pairs": {"title": "Both mates unmapped", "format": "{:,.0f}", "hidden": True},
    }
    content = {
        "id": "hgt_site",
        "section_name": "hgt-site",
        "description": "Candidate horizontal gene transfer sites called by hgt-site.",
        "plot_type": "table",
        "pconfig": {"id": "hgt_site_table", "title": "hgt-site: candidate HGT sites"},
        "headers": {k: v for k, v in headers.items() if k in row},
        "data": {sample: row},
    }
    with open(path, "w") as handle:
        json.dump(content, handle, indent=2)
    return row
