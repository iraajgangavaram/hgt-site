"""Measure how close hgt-site places each junction, with and without split reads.

Usage: python junctions.py <out.json> <run_dir> [<run_dir> ...]

Each run directory comes from run_benchmark.py --keep-sorted-bam and holds
host.sorted.bam, kraken2.txt and truth.bed.
"""

import json
import sys
import time
from pathlib import Path

from hgt_site.call import call_sites
from hgt_site.sites import site_filter

NODES = "/home/claude/bench/data/k2db/taxonomy/nodes.dmp"


def pair_estimate(st):
    for key in ("left_junction", "right_junction"):
        if st[key] != "":
            return int(st[key])
    return None


def split_estimate(st):
    for key in ("left_breakpoint", "right_breakpoint"):
        if st[key] != "":
            return int(st[key])
    return None


def evaluate(run):
    run = Path(run)
    truth = [(c, int(s)) for c, s, *_ in (line.split("\t") for line in open(run / "truth.bed"))]
    args = (str(run / "host.sorted.bam"), str(run / "kraken2.txt"), 9606, NODES)
    t0 = time.perf_counter()
    without = call_sites(*args, split_reads=False)
    t1 = time.perf_counter()
    sites = call_sites(*args)
    t2 = time.perf_counter()
    rows = []
    for site in sites:
        if site_filter(site) != "PASS":
            continue
        near = [p for c, p in truth if c == site.chrom and site.start - 500 <= p <= site.end + 500]
        if not near:
            continue
        true = min(near, key=lambda p: abs(p - site.start))
        st = site.stats
        split = split_estimate(st)
        rows.append({
            "truth": true,
            "pairs_error": abs(pair_estimate(st) - true),
            "split_error": abs(split - true) if split is not None else None,
            "split_reads": st["left_split_reads"] + st["right_split_reads"],
        })
    passed_with = sum(site_filter(s) == "PASS" for s in sites)
    passed_without = sum(site_filter(s) == "PASS" for s in without)
    return {
        "run": run.name,
        "sites": rows,
        "passed": passed_with,
        "passed_without_split_reads": passed_without,
        "seconds_with": round(t2 - t1, 2),
        "seconds_without": round(t1 - t0, 2),
    }


if __name__ == "__main__":
    out, *runs = sys.argv[1:]
    results = [evaluate(r) for r in runs]
    Path(out).write_text(json.dumps(results, indent=1))
    for r in results:
        errs = [s["split_error"] for s in r["sites"] if s["split_error"] is not None]
        exact = sum(e == 0 for e in errs)
        print(r["run"], len(r["sites"]), "matched;", len(errs), "with split reads;", exact, "exact;",
              "pairs median", sorted(s["pairs_error"] for s in r["sites"])[len(r["sites"]) // 2] if r["sites"] else "-",
              "| time", r["seconds_without"], "->", r["seconds_with"])
