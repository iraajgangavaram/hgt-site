"""Draw the README benchmark figures, in light and dark versions.

Usage: python figures.py <accuracy.json> <scaling.json> <junctions.json> <outdir>

accuracy.json:  list of {sim, cov, chim, seed, method, precision, recall, calls}
scaling.json:   list of {records, version, median}
junctions.json: output of junctions.py
"""

import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

THEMES = {
    "light": dict(text="#0b0b0b", muted="#52514e", grid="#e4e3df",
                  hgt="#2a78d6", lgt="#eb6834", sup="#1baf7a", alt="#eda100", tint="#a9c9ef"),
    "dark": dict(text="#f0efec", muted="#c3c2b7", grid="#33332f",
                 hgt="#3987e5", lgt="#d95926", sup="#199e70", alt="#c98500", tint="#2b4f7a"),
}
METHODS = [("hgt_site", "hgt-site", "hgt", "o"), ("lgtseek", "LGTSeek", "lgt", "s"),
           ("support_only", "Read support only", "sup", "^")]
CHIM = [2000, 10000, 30000]


def style(ax, c):
    ax.set_facecolor("none")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(c["grid"])
    ax.tick_params(colors=c["muted"], length=0, labelsize=9.5)
    ax.grid(axis="y", color=c["grid"], lw=0.8)
    ax.set_axisbelow(True)


def save(fig, out, name, mode):
    fig.patch.set_alpha(0)
    fig.savefig(Path(out) / f"{name}_{mode}.svg", transparent=True)
    plt.close(fig)


def accuracy(rows, out):
    by = defaultdict(list)
    for r in rows:
        by[(r["cov"], r["chim"], r["method"])].append(r)
    for mode, c in THEMES.items():
        fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.2), dpi=100, sharex=True)
        for col, cov in enumerate((10.0, 5.0)):
            for row, metric in enumerate(("precision", "recall")):
                ax = axes[row][col]
                style(ax, c)
                for offset, (key, label, colour, marker) in zip((-0.12, 0, 0.12), METHODS):
                    means = []
                    for x, chim in enumerate(CHIM):
                        vals = [r[metric] for r in by[(cov, chim, key)]]
                        ax.scatter([x + offset] * len(vals), vals, s=14, color=c[colour], alpha=0.35,
                                   marker=marker, linewidths=0, zorder=2)
                        means.append(mean(vals))
                    ax.plot([x + offset for x in range(3)], means, color=c[colour], lw=2, marker=marker,
                            ms=8, label=label, zorder=3)
                if metric == "precision":
                    ax.set_yscale("log")
                    ax.set_ylim(0.004, 1.6)
                    ax.set_yticks([0.01, 0.1, 1], ["1%", "10%", "100%"])
                else:
                    ax.set_ylim(0.7, 1.03)
                    ax.set_yticks([0.7, 0.8, 0.9, 1.0], ["70%", "80%", "90%", "100%"])
                ax.minorticks_off()
                if row == 0:
                    ax.set_title(f"{cov:g}x coverage", color=c["text"], fontsize=11, loc="left")
                if col == 0:
                    ax.set_ylabel("Precision (log scale)" if metric == "precision" else "Recall",
                                  color=c["muted"], fontsize=10)
                if row == 1:
                    ax.set_xticks(range(3), [f"{x:,}" for x in CHIM])
                    ax.set_xlabel("Chimeric read pairs in the sample", color=c["muted"], fontsize=10)
        leg = axes[0][1].legend(loc="lower left", frameon=False, fontsize=9.5)
        for text in leg.get_texts():
            text.set_color(c["text"])
        fig.tight_layout(h_pad=1.5)
        save(fig, out, "benchmark_accuracy", mode)


def speed(out):
    tools = [("hgt-site", 12.0, "Kraken2 on unmapped reads", 2.7, "hgt"),
             ("LGTSeek", 129.0, "Second alignment to bacterial genomes", 25.0, "lgt")]
    for mode, c in THEMES.items():
        fig, ax = plt.subplots(figsize=(8.2, 2.6), dpi=100)
        style(ax, c)
        ax.grid(axis="x", color=c["grid"], lw=0.8)
        ax.grid(axis="y", visible=False)
        for y, (_name, prep, prep_label, call, colour) in enumerate(tools):
            ax.barh(y, prep, color=c[colour], alpha=0.35, height=0.42, edgecolor="none")
            ax.barh(y, call, left=prep + 0.6, color=c[colour], height=0.42, edgecolor="none")
            ax.text(prep + call + 3, y, f"{prep + call:.0f} s", va="center", color=c["text"], fontsize=10.5,
                    fontweight="bold")
            ax.text(0.5, y + 0.3, f"{prep_label} {prep:.0f} s, then calling {call:g} s",
                    va="top", color=c["muted"], fontsize=9)
        ax.set_yticks([0, 1], [t[0] for t in tools], color=c["text"], fontsize=10.5)
        ax.set_ylim(-0.4, 1.75)
        ax.invert_yaxis()
        ax.set_xlim(0, 175)
        ax.set_xlabel("Seconds after the host alignment, 2.1 million reads, 2 CPU cores", color=c["muted"],
                      fontsize=10)
        fig.tight_layout()
        save(fig, out, "benchmark_speed", mode)


def scaling(points, out):
    series = defaultdict(list)
    for p in points:
        series[p["version"]].append((p["records"] / 1e6, p["median"]))
    for mode, c in THEMES.items():
        fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=100)
        style(ax, c)
        for version, colour, label in (("0.1.1", "alt", "version 0.1.1"), ("current", "hgt", "current version")):
            xy = sorted(series[version])
            ax.plot([x for x, _ in xy], [y for _, y in xy], color=c[colour], lw=2, marker="o", ms=8)
            x, y = xy[-1]
            ax.text(x + 0.25, y, label, color=c["text"], va="center", fontsize=10, fontweight="bold")
        ax.set_xlim(0, max(x for x, _ in series["current"]) * 1.32)
        ax.set_ylim(0, None)
        ax.set_xlabel("Records in the BAM file (millions)", color=c["muted"], fontsize=10)
        ax.set_ylabel("Seconds (median of 3 runs)", color=c["muted"], fontsize=10)
        ax.set_title("hgt-site run time grows linearly with file size", color=c["text"], fontsize=11, loc="left")
        fig.tight_layout()
        save(fig, out, "benchmark_scaling", mode)


def fragments(out):
    sizes = ["400 bp", "800 bp", "1000 bp"]
    fixed = [100 / 100, 201 / 100, 239 / 100]
    auto = [100 / 100, 113 / 100, 121 / 100]
    for mode, c in THEMES.items():
        fig, ax = plt.subplots(figsize=(7.2, 3.4), dpi=100)
        style(ax, c)
        w = 0.36
        for i, (f, a) in enumerate(zip(fixed, auto)):
            ax.bar(i - w / 2 - 0.01, f, w, color=c["alt"], edgecolor="none")
            ax.bar(i + w / 2 + 0.01, a, w, color=c["hgt"], edgecolor="none")
            ax.text(i - w / 2 - 0.01, f + 0.05, f"{f:.2f}", ha="center", color=c["text"], fontsize=9.5)
            ax.text(i + w / 2 + 0.01, a + 0.05, f"{a:.2f}", ha="center", color=c["text"], fontsize=9.5)
        ideal = ax.axhline(1, color=c["muted"], lw=1, ls=(0, (4, 3)))
        ax.set_xticks(range(3), sizes)
        ax.set_ylim(0, 2.75)
        ax.set_ylabel("Calls per real insertion", color=c["muted"], fontsize=10)
        ax.set_xlabel("Mean fragment length of the library", color=c["muted"], fontsize=10)
        handles = [plt.Rectangle((0, 0), 1, 1, color=c["alt"]), plt.Rectangle((0, 0), 1, 1, color=c["hgt"])]
        leg = ax.legend(handles + [ideal], ["Fixed clustering limits", "Limits from measured fragment size",
                                            "Ideal: one call per insertion"],
                        loc="upper left", frameon=False, fontsize=9.5)
        for text in leg.get_texts():
            text.set_color(c["text"])
        fig.tight_layout()
        save(fig, out, "benchmark_fragments", mode)


def junctions(runs, out):
    """How close each junction is, with read pairs alone and with split reads."""
    sites = [s for run in runs for s in run["sites"]]
    pairs = sorted(s["pairs_error"] for s in sites)
    split = sorted(s["split_error"] if s["split_error"] is not None else s["pairs_error"] for s in sites)
    n = len(sites)
    cuts = [0, 1, 2, 5, 10, 25, 50]

    def share(values):
        return [100 * sum(v <= c for v in values) / n for c in cuts]

    for mode, c in THEMES.items():
        fig, ax = plt.subplots(figsize=(7.2, 3.6), dpi=100)
        style(ax, c)
        x = range(len(cuts))
        w = 0.36
        for offset, values, colour, label in (
            (-w / 2 - 0.01, pairs, "alt", "Read pairs only"),
            (w / 2 + 0.01, split, "hgt", "With split reads"),
        ):
            heights = share(values)
            ax.bar([i + offset for i in x], heights, w, color=c[colour], edgecolor="none", label=label)
            for i, h in enumerate(heights):
                ax.text(i + offset, h + 1.5, f"{h:.0f}", ha="center", color=c["text"], fontsize=9)
        ax.set_xticks(list(x), ["exact"] + [f"{v} bp" for v in cuts[1:]])
        ax.set_ylim(0, 126)
        ax.set_yticks([0, 25, 50, 75, 100], ["0", "25%", "50%", "75%", "100%"])
        ax.set_ylabel("Insertions placed this close", color=c["muted"], fontsize=10)
        ax.set_xlabel("Distance from the true junction", color=c["muted"], fontsize=10)
        ax.set_title(f"Junction accuracy, {n} insertions", color=c["text"], fontsize=11, loc="left")
        leg = ax.legend(loc="upper right", bbox_to_anchor=(1, 1.04), ncol=2, frameon=False, fontsize=9.5)
        for text in leg.get_texts():
            text.set_color(c["text"])
        fig.tight_layout()
        save(fig, out, "benchmark_junctions", mode)


if __name__ == "__main__":
    acc, scale, junc, outdir = sys.argv[1:5]
    Path(outdir).mkdir(parents=True, exist_ok=True)
    accuracy(json.load(open(acc)), outdir)
    speed(outdir)
    scaling(json.load(open(scale)), outdir)
    fragments(outdir)
    junctions(json.load(open(junc)), outdir)
