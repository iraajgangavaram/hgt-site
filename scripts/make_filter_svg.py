"""Build the animated figure of what the enrichment test removes.

Three stages, from one real benchmark run: every position where foreign reads
cluster, then the subset with enough read support to look convincing, then the
few that survive the test against chance chimeras.

The run used is the hardest of the benchmark grid, 30,000 chance chimeras at
ten-fold coverage, because that is the condition where read support genuinely
fails: six rejected candidates carry as much support as the weakest real
insertion, so no threshold on read count separates them. On an easy run the
figure would be dishonest, since counting reads would almost do the job on its
own.

The animation carries the argument. A still of any one stage looks
unremarkable; it is the collapse from one to the next that shows what the tool
does.

    python scripts/make_filter_svg.py path/to/sites.tsv

Defaults to the run committed under benchmark/figure_data/, gzipped or not.
Every candidate with three or more supporting reads is drawn. The cloud below
that is drawn as a random sample, because eight thousand overlapping dots cost
a quarter of a megabyte and show nothing the sample does not; the counts in the
captions are the real ones, and the caption in the README says the cloud is
sampled.

Animation is plain CSS inside the SVG, because GitHub renders a README image
in an <img> tag, where scripts never run but declarative animation does.
"""

import csv
import gzip
import math
import random
import sys
from pathlib import Path

WIDTH, HEIGHT = 920, 340
LEFT, RIGHT, TOP, BOTTOM = 64, 26, 74, 54
SUPPORT_LINE = 3
CLOUD_SAMPLE = 1800
LOOP = 15.0

THEMES = {
    "light": {
        "axis": "#8d8a84", "label": "#272524", "faint": "#cac6bf", "grey": "#b9b5ae",
        "mid": "#c25c10", "pass": "#157954", "bg": "none", "note": "#6d6a65",
    },
    "dark": {
        "axis": "#7e8a97", "label": "#e6e4e0", "faint": "#3f444b", "grey": "#5a6068",
        "mid": "#f0a35e", "pass": "#4fd1a0", "bg": "none", "note": "#9b968e",
    },
}


def load(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [
        (
            (int(r["start"]) + int(r["end"])) // 2,
            int(r["n_pairs"]),
            r["filter"] == "PASS",
        )
        for r in rows
    ]


def build(points, theme_name):
    t = THEMES[theme_name]
    rng = random.Random(3)
    span = max(p[0] for p in points)
    top_pairs = max(p[1] for p in points)
    plot_w = WIDTH - LEFT - RIGHT
    plot_h = HEIGHT - TOP - BOTTOM

    def x_of(pos):
        return LEFT + plot_w * pos / span

    def y_of(pairs):
        # Square root keeps the crowded low counts legible without pretending
        # the axis is linear.
        return HEIGHT - BOTTOM - plot_h * math.sqrt(pairs) / math.sqrt(top_pairs)

    cloud = [p for p in points if p[1] < SUPPORT_LINE]
    kept = [p for p in points if p[1] >= SUPPORT_LINE]
    hits = [p for p in points if p[2]]
    rejected = [p for p in kept if not p[2]]
    # The honest point of the whole figure: read support does not separate
    # these, because some rejected candidates carry as much of it as the
    # weakest real site.
    weakest = min(p[1] for p in hits)
    as_well_supported = sum(1 for p in rejected if p[1] >= weakest)

    shown = cloud if len(cloud) <= CLOUD_SAMPLE else rng.sample(cloud, CLOUD_SAMPLE)

    def dots(group, radius):
        return " ".join(
            f"M{x_of(pos) + rng.uniform(-0.7, 0.7):.1f} "
            f"{y_of(pairs) + rng.uniform(-1.5, 1.5):.1f}h{radius:.1f}"
            for pos, pairs, _ in group
        )

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" '
        f"font-family=\"-apple-system,'Segoe UI',Helvetica,Arial,sans-serif\" role=\"img\" "
        f'aria-label="A benchmark run with 30,000 chance chimeras gives {len(points):,} '
        f'positions where foreign reads cluster. {len(kept):,} of them have {SUPPORT_LINE} '
        f'or more supporting reads, and {as_well_supported} of those have as much support as '
        f'the weakest real insertion, so read count alone cannot separate them. '
        f'{len(hits)} survive the test against chance chimeras and every one is a real '
        f'insertion. The low-support cloud is drawn as a random sample.">',
        "<style>",
        # Stage 1: everything, grey. Stage 2: the thinly supported go, and what
        # is left turns orange. Stage 3: all but the survivors go.
        "@keyframes cloud{0%{opacity:0}5%,32%{opacity:.9}40%,100%{opacity:0}}",
        f"@keyframes kept{{0%{{opacity:0;stroke:{t['grey']}}}"
        f"5%,33%{{opacity:1;stroke:{t['grey']}}}"
        f"40%,68%{{opacity:1;stroke:{t['mid']}}}"
        f"76%,100%{{opacity:0;stroke:{t['mid']}}}}}",
        "@keyframes hits{0%,74%{opacity:0;r:3}81%,100%{opacity:1;r:5}}",
        "@keyframes ghost{0%,76%{opacity:0}84%,100%{opacity:1}}",
        "@keyframes cap1{0%{opacity:0}5%,31%{opacity:1}37%,100%{opacity:0}}",
        "@keyframes cap2{0%,37%{opacity:0}43%,68%{opacity:1}74%,100%{opacity:0}}",
        "@keyframes cap3{0%,74%{opacity:0}81%,100%{opacity:1}}",
    ]
    for name in ("cloud", "kept", "hits", "ghost", "cap1", "cap2", "cap3"):
        ease = "ease-out" if name == "hits" else "ease-in-out"
        out.append(f".{name}{{animation:{name} {LOOP}s {ease} infinite}}")
    out.append("</style>")

    # Axes.
    out.append(
        f'<line x1="{LEFT}" y1="{HEIGHT - BOTTOM}" x2="{WIDTH - RIGHT}" y2="{HEIGHT - BOTTOM}" '
        f'stroke="{t["axis"]}" stroke-width="1.2"/>'
    )
    out.append(
        f'<line x1="{LEFT}" y1="{TOP - 10}" x2="{LEFT}" y2="{HEIGHT - BOTTOM}" '
        f'stroke="{t["axis"]}" stroke-width="1.2"/>'
    )
    for pairs in (1, 5, 10, 20, 50):
        if pairs > top_pairs:
            continue
        y = y_of(pairs)
        out.append(
            f'<text x="{LEFT - 10}" y="{y + 4:.1f}" font-size="11" fill="{t["note"]}" '
            f'text-anchor="end">{pairs}</text>'
        )
    out.append(
        f'<text x="15" y="{TOP + plot_h / 2:.0f}" font-size="11.5" fill="{t["note"]}" '
        f'text-anchor="middle" '
        f'transform="rotate(-90 15 {TOP + plot_h / 2:.0f})">supporting reads</text>'
    )
    out.append(
        f'<text x="{LEFT + plot_w / 2:.0f}" y="{HEIGHT - 16}" font-size="11.5" '
        f'fill="{t["note"]}" text-anchor="middle">position along 10 Mb of host genome</text>'
    )

    out.append(
        f'<path class="cloud" d="{dots(shown, 1.6)}" stroke="{t["faint"]}" '
        f'stroke-width="2.3" stroke-linecap="round" fill="none"/>'
    )
    out.append(
        f'<path class="kept" d="{dots(kept, 1.8)}" stroke="{t["grey"]}" '
        f'stroke-width="2.7" stroke-linecap="round" fill="none"/>'
    )
    # The rejected candidates that carry as much support as the weakest real
    # site stay on in the last stage, as rings. Stating in a caption that read
    # count cannot separate them is weaker than showing where they sit.
    for pos, pairs, _ in rejected:
        if pairs < weakest:
            continue
        out.append(
            f'<circle class="ghost" cx="{x_of(pos):.1f}" cy="{y_of(pairs):.1f}" r="5" '
            f'fill="none" stroke="{t["grey"]}" stroke-width="1.6"/>'
        )
    for pos, pairs, _ in hits:
        out.append(
            f'<circle class="hits" cx="{x_of(pos):.1f}" cy="{y_of(pairs):.1f}" r="3" '
            f'fill="{t["pass"]}"/>'
        )
    # Low and to the right, where nothing is left by the last stage.
    legend_y = HEIGHT - BOTTOM - 34
    out.append(
        f'<g class="ghost">'
        f'<circle cx="{WIDTH - RIGHT - 232}" cy="{legend_y}" r="5" fill="{t["pass"]}"/>'
        f'<text x="{WIDTH - RIGHT - 220}" y="{legend_y + 4}" font-size="11.5" '
        f'fill="{t["note"]}">real insertion</text>'
        f'<circle cx="{WIDTH - RIGHT - 128}" cy="{legend_y}" r="5" fill="none" '
        f'stroke="{t["grey"]}" stroke-width="1.6"/>'
        f'<text x="{WIDTH - RIGHT - 116}" y="{legend_y + 4}" font-size="11.5" '
        f'fill="{t["note"]}">chance chimera, rejected</text>'
        f"</g>"
    )

    captions = [
        ("cap1", f"{len(points):,} positions where foreign reads cluster", None),
        (
            "cap2",
            f"{len(kept):,} have {SUPPORT_LINE} or more supporting reads",
            f"{as_well_supported} of them carry as much support as the weakest real "
            f"insertion, so counting reads cannot separate them",
        ),
        (
            "cap3",
            f"{len(hits)} survive the test against chance chimeras",
            "every one of them is a real insertion",
        ),
    ]
    for cls, headline, note in captions:
        out.append(
            f'<g class="{cls}"><text x="{LEFT}" y="{TOP - 38}" font-size="15.5" '
            f'font-weight="600" fill="{t["label"]}">{headline}</text>'
        )
        if note:
            out.append(
                f'<text x="{LEFT}" y="{TOP - 18}" font-size="12.5" '
                f'fill="{t["note"]}">{note}</text>'
            )
        out.append("</g>")
    out.append("</svg>")
    return "\n".join(out)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    default = root / "benchmark" / "figure_data" / "sites.tsv.gz"
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else default
    points = load(source)
    for name in THEMES:
        path = root / "docs" / "images" / f"filter-{name}.svg"
        path.write_text(build(points, name))
        print(f"wrote {path} ({path.stat().st_size:,} bytes)")
