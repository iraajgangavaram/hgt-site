"""Compare called sites with a truth BED file."""

import csv


def read_truth(path):
    sites = []
    with open(path) as handle:
        for line in handle:
            cols = line.split("\t")
            sites.append((cols[0], int(cols[1])))
    return sites


def read_calls(path, passed_only=True):
    calls = []
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if passed_only and row["filter"] != "PASS":
                continue
            calls.append((row["chrom"], int(row["start"]), int(row["end"])))
    return calls


def evaluate(calls, truth, slop=500):
    """Match calls to true sites; a call matches if the site lies within slop bp of it."""
    matched_truth = set()
    true_calls = 0
    for chrom, start, end in calls:
        hits = [
            i
            for i, (tchrom, tpos) in enumerate(truth)
            if tchrom == chrom and start - slop <= tpos <= end + slop
        ]
        if hits:
            true_calls += 1
            matched_truth.update(hits)
    precision = true_calls / len(calls) if calls else 0.0
    recall = len(matched_truth) / len(truth) if truth else 0.0
    return {
        "calls": len(calls),
        "true_calls": true_calls,
        "truth": len(truth),
        "recovered": len(matched_truth),
        "precision": precision,
        "recall": recall,
    }
