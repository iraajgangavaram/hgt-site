"""Combine the sites called in several samples into one table.

Calling samples one at a time leaves the most useful comparison undone. A
position that turns up in every sample, including controls, is almost always a
property of the protocol, the reference or the reagents rather than of the
biology, however convincing its statistics look in any one sample. A position
in one sample of a group, with split reads behind it, is the opposite.

This module takes the per-sample tables hgt-site already writes, merges
positions that fall close together, and reports how many samples each one
appeared in. It reads the output of calling rather than the alignments, so a
cohort can be assembled from runs made at different times, and recalling is
never needed.
"""

import csv
import logging
from pathlib import Path

log = logging.getLogger("hgt_site")

# Positions this far apart or closer are treated as the same site across
# samples. Junctions estimated from read pairs can sit a fragment away from
# the truth, so a few hundred bases is the right scale.
WINDOW = 500

# A site seen in at least this share of the samples is flagged as common,
# which usually means an artefact rather than a shared insertion.
COMMON_FRACTION = 0.8

COLUMNS = [
    "chrom",
    "start",
    "end",
    "n_samples",
    "n_samples_passed",
    "samples",
    "samples_passed",
    "top_taxid",
    "total_pairs",
    "max_pairs",
    "best_qvalue",
    "split_reads",
    "gene",
    "flags",
]


def read_sites(path, include_failed=False):
    """Read one hgt-site table. Returns the rows, as dicts, in file order."""
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "chrom" not in reader.fieldnames:
            raise ValueError(f"{path} does not look like an hgt-site table: no chrom column.")
        for row in reader:
            if not include_failed and row.get("filter") != "PASS":
                continue
            rows.append(row)
    return rows


def _number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _sample_name(path):
    """A readable name for a sample, from its file name."""
    name = Path(path).name
    for suffix in (".tsv", ".txt", ".hgt-site", ".sites"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name or str(path)


def merge(samples, window=WINDOW, common_fraction=COMMON_FRACTION):
    """Merge the sites of several samples into one list of rows.

    samples maps a sample name to its rows. Sites are merged when they sit on
    the same contig within window bases of one another, which is deliberately
    taxon-blind: two samples reporting different taxa at the same position is
    itself a sign of a mismapping-prone region, and hiding that by keeping
    them apart would be unhelpful.
    """
    total = len(samples)
    entries = []
    for name, rows in samples.items():
        for row in rows:
            entries.append((row["chrom"], int(row["start"]), int(row["end"]), name, row))
    entries.sort(key=lambda e: (e[0], e[1], e[2]))

    merged = []
    for chrom, start, end, name, row in entries:
        if merged and merged[-1]["chrom"] == chrom and start - merged[-1]["end"] <= window:
            cluster = merged[-1]
            cluster["end"] = max(cluster["end"], end)
        else:
            cluster = {"chrom": chrom, "start": start, "end": end, "members": []}
            merged.append(cluster)
        cluster["members"].append((name, row))

    out = []
    for cluster in merged:
        members = cluster["members"]
        names = sorted({name for name, _ in members})
        passed = sorted({name for name, row in members if row.get("filter") == "PASS"})
        taxa = {}
        for _, row in members:
            taxa[row.get("top_taxid", "")] = taxa.get(row.get("top_taxid", ""), 0) + 1
        pairs = [int(_number(row.get("n_pairs"))) for _, row in members]
        qvalues = [_number(row.get("qvalue"), 1.0) for _, row in members]
        splits = sum(
            int(_number(row.get("left_split_reads"))) + int(_number(row.get("right_split_reads")))
            for _, row in members
        )
        genes = sorted({row.get("gene", "") for _, row in members if row.get("gene")})

        flags = []
        if total > 1 and len(names) >= max(2, round(common_fraction * total)):
            flags.append("COMMON")
        if len(names) > 1:
            flags.append("RECURRENT")
        if len(taxa) > 1:
            flags.append("TAXA_DIFFER")

        out.append({
            "chrom": cluster["chrom"],
            "start": cluster["start"],
            "end": cluster["end"],
            "n_samples": len(names),
            "n_samples_passed": len(passed),
            "samples": ";".join(names),
            "samples_passed": ";".join(passed) or ".",
            "top_taxid": max(taxa, key=lambda t: taxa[t]) if taxa else "",
            "total_pairs": sum(pairs),
            "max_pairs": max(pairs) if pairs else 0,
            "best_qvalue": f"{min(qvalues):.4g}" if qvalues else "",
            "split_reads": splits,
            "gene": ";".join(genes),
            "flags": ";".join(flags) if flags else ".",
        })
    out.sort(key=lambda r: (r["chrom"], r["start"]))
    return out


def write_cohort(rows, path):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return len(rows)


def write_matrix(rows, samples, path):
    """A position by sample grid of supporting pair counts, for a quick look."""
    names = sorted(samples)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["site", *names])
        for row in rows:
            site = f"{row['chrom']}:{row['start']}-{row['end']}"
            present = set(row["samples"].split(";"))
            writer.writerow([site, *("1" if name in present else "0" for name in names)])
    return len(rows)


def build(paths, names=None, window=WINDOW, common_fraction=COMMON_FRACTION, include_failed=False):
    """Read per-sample tables and merge them. Returns (rows, sample names)."""
    if names and len(names) != len(paths):
        raise ValueError(
            f"{len(names)} sample names were given for {len(paths)} files. "
            "Give one name per file, in the same order, or none at all."
        )
    samples = {}
    for index, path in enumerate(paths):
        name = names[index] if names else _sample_name(path)
        if name in samples:
            raise ValueError(f"Two samples are both called {name!r}. Use --names to tell them apart.")
        samples[name] = read_sites(path, include_failed)
        log.info("%s: %d sites.", name, len(samples[name]))
    rows = merge(samples, window, common_fraction)
    log.info(
        "%d distinct positions across %d samples; %d seen in more than one.",
        len(rows), len(samples), sum(1 for r in rows if r["n_samples"] > 1),
    )
    return rows, list(samples)
