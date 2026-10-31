"""Pool benchmark summary.json files into one table per setting."""

import json
import sys
from collections import defaultdict
from pathlib import Path


def main(runs_dir):
    groups = defaultdict(list)
    for path in sorted(Path(runs_dir).glob("*/summary.json")):
        s = json.loads(path.read_text())
        groups[(s["coverage"], s["chimera_pairs"])].append(s)

    methods = ["hgt_site", "support_only"]
    if all("lgtseek" in r for runs in groups.values() for r in runs):
        methods.append("lgtseek")
    labels = {"hgt_site": "hgt-site", "support_only": "Support-only", "lgtseek": "LGTSeek"}
    header = ["Coverage", "Chimeric pairs", "Runs"]
    for m in methods:
        header += [f"{labels[m]} precision", f"{labels[m]} recall"]
    print("| " + " | ".join(header) + " |")
    print("|" + "---|" * len(header))
    for (cov, chim), runs in sorted(groups.items(), key=lambda kv: (-kv[0][0], kv[0][1])):
        row = [f"{cov:g}x", f"{chim:,}", str(len(runs))]
        for method in methods:
            calls = sum(r[method]["calls"] for r in runs)
            true_calls = sum(r[method]["true_calls"] for r in runs)
            truth = sum(r[method]["truth"] for r in runs)
            recovered = sum(r[method]["recovered"] for r in runs)
            row += [
                f"{true_calls / calls:.3f} ({true_calls}/{calls})" if calls else "n/a",
                f"{recovered / truth:.3f} ({recovered}/{truth})",
            ]
        print("| " + " | ".join(row) + " |")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "runs")
