"""Draw a terminal session as an SVG, in light and dark versions.

Hand-drawn rather than recorded, so the text is real output taken from a run
and the picture stays small and readable in a README.
"""

from pathlib import Path
from xml.sax.saxutils import escape

OUT = Path("/home/claude/hgt-site-repo/docs/images")

THEMES = {
    "light": dict(bg="#faf9f7", frame="#e4e3df", text="#0b0b0b", muted="#52514e",
                  prompt="#1baf7a", flag="#2a78d6", warn="#eda100", title="#52514e"),
    "dark": dict(bg="#1b1b19", frame="#33332f", text="#f0efec", muted="#c3c2b7",
                 prompt="#199e70", flag="#3987e5", warn="#c98500", title="#c3c2b7"),
}

# (style, text). Styles: prompt, out, muted, flag
LINES = [
    ("prompt", "pip install hgt-site"),
    ("prompt", "hgt-site run --bam sample.bam --kraken-db kraken2-db \\"),
    ("cont", "    --host-taxid 9606 --threads 8 -o sites.tsv --vcf sites.vcf"),
    ("blank", ""),
    ("muted", "INFO: Wrote 105187 unmapped reads."),
    ("muted", "INFO: Running Kraken2 on the unmapped reads."),
    ("muted", "INFO: 2107 anchored pairs have a Kraken2 call; 2107 are from taxa outside the host."),
    ("muted", "INFO: Fragment size from 100000 pairs: median 400 bp, 99th percentile 494 bp."),
    ("muted", "INFO: Exact junctions from split reads found for 20 of 21 sites that passed."),
    ("flag", "1701 candidate sites, 21 passed filters -> sites.tsv"),
    ("blank", ""),
    ("prompt", "head -3 sites.tsv | cut -f1-3,6,7,10,13"),
    ("head", "chrom        start     end       left_breakpoint  split_reads  n_pairs  top_taxid"),
    ("out", "chr21_20M    638838    639307    639072           8            11       562"),
    ("out", "chr21_20M    819852    820366    820112           4            18       562"),
]

CHAR_W, LINE_H, PAD, TOP = 7.8, 21.0, 20.0, 46.0
WIDTH = 860


def draw(mode):
    c = THEMES[mode]
    height = TOP + LINE_H * len(LINES) + PAD
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{height:.0f}" '
        f'viewBox="0 0 {WIDTH} {height:.0f}" role="img" aria-label="A terminal showing hgt-site being '
        f'installed and run, printing 21 candidate sites that passed the filters">',
        f'<rect x="0.5" y="0.5" width="{WIDTH - 1}" height="{height - 1:.0f}" rx="10" fill="{c["bg"]}" '
        f'stroke="{c["frame"]}"/>',
        f'<line x1="0" y1="32" x2="{WIDTH}" y2="32" stroke="{c["frame"]}"/>',
    ]
    for i, colour in enumerate(("#eb6834", "#eda100", "#1baf7a")):
        parts.append(f'<circle cx="{22 + i * 18}" cy="16.5" r="5" fill="{colour}" opacity="0.85"/>')
    parts.append(
        f'<text x="{WIDTH / 2}" y="21" text-anchor="middle" font-family="Inter, system-ui, sans-serif" '
        f'font-size="11.5" fill="{c["title"]}">hgt-site</text>'
    )
    font = 'font-family="ui-monospace, SFMono-Regular, Menlo, Consolas, monospace" font-size="12.5"'
    y = TOP + 14
    for style, text in LINES:
        if style == "blank":
            y += LINE_H
            continue
        x = PAD
        if style == "prompt":
            parts.append(f'<text x="{x}" y="{y:.0f}" {font} fill="{c["prompt"]}">$</text>')
            x += 2 * CHAR_W
            fill, weight = c["text"], ' font-weight="600"'
        elif style == "cont":
            x += 2 * CHAR_W
            fill, weight = c["text"], ' font-weight="600"'
        elif style == "muted":
            fill, weight = c["muted"], ""
        elif style == "head":
            fill, weight = c["muted"], ""
        elif style == "flag":
            fill, weight = c["flag"], ' font-weight="600"'
        else:
            fill, weight = c["text"], ""
        parts.append(
            f'<text x="{x:.0f}" y="{y:.0f}" {font}{weight} fill="{fill}" '
            f'xml:space="preserve">{escape(text)}</text>'
        )
        y += LINE_H
    parts.append("</svg>")
    path = OUT / f"demo_{mode}.svg"
    path.write_text("\n".join(parts) + "\n")
    return path


for mode in THEMES:
    print(draw(mode))
