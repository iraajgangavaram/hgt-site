"""Build the animated terminal demo shown in the README.

The session is a real one: the commands and the output below were produced by
running hgt-site on its own simulated demo dataset with seed 1. Regenerate the
pair with:

    python scripts/make_demo_svg.py

Animation is plain CSS inside the SVG, because GitHub renders a README image
in an <img> tag, where scripts never run but declarative animation does.
"""

from pathlib import Path

# Each line is (delay in seconds, kind, text). "cmd" lines are typed out, the
# rest appear at once, as a terminal does.
LINES = [
    (0.0, "cmd", "conda install -c conda-forge -c bioconda hgt-site"),
    (1.8, "dim", "done"),
    (2.4, "cmd", "hgt-site simulate -o demo --seed 1"),
    (4.0, "dim", "Wrote demo.sam, .kraken2.txt, .truth.bed and .host.fa (5 insertions)"),
    (4.8, "cmd", "hgt-site call --bam demo.sam --kraken demo.kraken2.txt \\"),
    (6.2, "cmd2", "              --host-taxid 9606 -o sites.tsv"),
    (7.2, "dim", "281 anchored pairs, 234 with a taxonomic call outside the host"),
    (7.9, "warn", "112 candidate sites, 5 passed filters -> sites.tsv"),
    (8.8, "head", "chrom      start     end  pairs   taxid  filter"),
    (9.2, "row", "chr1       31897   32626     20     562  PASS"),
    (9.5, "row", "chr1       48529   49096     24     562  PASS"),
    (9.8, "row", "chr1       93829   94348     13     562  PASS"),
    (10.1, "row", "chr2       30954   31468     23     562  PASS"),
    (10.4, "row", "chr2       32296   32878     21     562  PASS"),
    (11.4, "note", "107 rejected, each with the reason in the filter column"),
]

LOOP = 16.0
CHAR = 8.4
LINE_HEIGHT = 25
TOP = 74
PAD = 26
WIDTH = 920

THEMES = {
    "light": {
        "bg": "#ffffff", "chrome": "#f2f1ef", "border": "#dcdad6", "text": "#2b2a28",
        "dim": "#6d6a65", "prompt": "#1b7f5a", "cmd": "#2b2a28", "warn": "#b4540a",
        "pass": "#1b7f5a", "head": "#4a4845", "shadow": "0.10",
    },
    "dark": {
        "bg": "#14161a", "chrome": "#1d2026", "border": "#2b2f36", "text": "#e3e1dd",
        "dim": "#9b968e", "prompt": "#4fd1a0", "cmd": "#e3e1dd", "warn": "#f0a35e",
        "pass": "#4fd1a0", "head": "#b9b5ae", "shadow": "0.45",
    },
}


def escape(text):
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def build(theme_name):
    t = THEMES[theme_name]
    height = TOP + LINE_HEIGHT * len(LINES) + PAD
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height}" '
        f'viewBox="0 0 {WIDTH} {height}" font-family="ui-monospace,SFMono-Regular,Menlo,Consolas,monospace" '
        f'font-size="14.5" role="img" aria-label="hgt-site terminal session: install, simulate a dataset, '
        f'call sites, and show five passing sites out of 112 candidates">',
        "<style>",
        ".ln{opacity:0;animation:show %ss steps(1,end) infinite}" % LOOP,
        "@keyframes show{0%{opacity:0}100%{opacity:1}}",
        ".cur{animation:blink 1.06s steps(1,end) infinite}",
        "@keyframes blink{0%,50%{opacity:1}51%,100%{opacity:0}}",
    ]
    # One reveal keyframe set per line, so each appears at its own moment and
    # they all clear together when the loop restarts.
    for i, (delay, _, _) in enumerate(LINES):
        start = delay / LOOP * 100
        out.append(
            f"@keyframes r{i}{{0%,{start:.3f}%{{opacity:0}}{start + 0.001:.3f}%,100%{{opacity:1}}}}"
        )
        out.append(f".l{i}{{animation-name:r{i}}}")
    # Commands type themselves out by widening a clip rectangle.
    for i, (delay, kind, text) in enumerate(LINES):
        if kind not in ("cmd", "cmd2"):
            continue
        typing = max(len(text) * 0.022, 0.25)
        start = delay / LOOP * 100
        end = (delay + typing) / LOOP * 100
        out.append(
            f"@keyframes t{i}{{0%,{start:.3f}%{{width:0}}{end:.3f}%,100%{{width:{WIDTH}px}}}}"
        )
        out.append(f".c{i}{{animation:t{i} {LOOP}s linear infinite}}")
    out.append("</style>")
    out.append(
        f'<rect x="2" y="2" width="{WIDTH - 4}" height="{height - 4}" rx="11" fill="{t["bg"]}" '
        f'stroke="{t["border"]}" stroke-width="1.5"/>'
    )
    out.append(
        f'<path d="M2 13a11 11 0 0 1 11-11h{WIDTH - 26}a11 11 0 0 1 11 11v26H2z" fill="{t["chrome"]}"/>'
    )
    out.append(f'<line x1="2" y1="39" x2="{WIDTH - 2}" y2="39" stroke="{t["border"]}" stroke-width="1.5"/>')
    for i, colour in enumerate(("#ec6a5e", "#f4bf4f", "#61c554")):
        out.append(f'<circle cx="{26 + i * 20}" cy="21" r="6" fill="{colour}"/>')
    out.append(
        f'<text x="{WIDTH / 2}" y="26" fill="{t["dim"]}" font-size="12.5" text-anchor="middle">hgt-site</text>'
    )

    for i, (delay, kind, text) in enumerate(LINES):
        y = TOP + i * LINE_HEIGHT
        if kind in ("cmd", "cmd2"):
            prompt = "$ " if kind == "cmd" else "  "
            out.append(f'<g class="ln l{i}">')
            if kind == "cmd":
                out.append(f'<text x="{PAD}" y="{y}" fill="{t["prompt"]}">$</text>')
            out.append(f'<clipPath id="{theme_name}{i}"><rect class="c{i}" x="{PAD + len(prompt) * CHAR}" '
                       f'y="{y - 14}" height="20" width="0"/></clipPath>')
            out.append(
                f'<text xml:space="preserve" x="{PAD + len(prompt) * CHAR}" y="{y}" fill="{t["cmd"]}" '
                f'clip-path="url(#{theme_name}{i})">{escape(text)}</text>'
            )
            out.append("</g>")
        else:
            colour = {
                "dim": t["dim"], "warn": t["warn"], "row": t["pass"],
                "head": t["head"], "note": t["dim"],
            }[kind]
            weight = ' font-weight="600"' if kind in ("warn", "head") else ""
            style = ' font-style="italic"' if kind == "note" else ""
            out.append(
                f'<text xml:space="preserve" class="ln l{i}" x="{PAD}" y="{y}" fill="{colour}"{weight}{style}>'
                f"{escape(text)}</text>"
            )

    last_y = TOP + (len(LINES) - 1) * LINE_HEIGHT + LINE_HEIGHT
    out.append(f'<text x="{PAD}" y="{last_y}" fill="{t["prompt"]}">$</text>')
    out.append(
        f'<rect class="cur" x="{PAD + 2 * CHAR}" y="{last_y - 13}" width="9" height="17" fill="{t["dim"]}"/>'
    )
    out.append("</svg>")
    return "\n".join(out)


if __name__ == "__main__":
    images = Path(__file__).resolve().parent.parent / "docs" / "images"
    images.mkdir(parents=True, exist_ok=True)
    for name in THEMES:
        path = images / f"demo-{name}.svg"
        path.write_text(build(name))
        print(f"wrote {path} ({path.stat().st_size:,} bytes)")
