"""Animate the existing logo: the pieces arrive, then the name assembles.

Same design, colours and wording as docs/images/logo_light.svg and its dark
twin, which stay in the repository unchanged. The host segments slide in from
the sides, the inserted segment drops in between them, the junction markers
draw themselves, and the letters of the name fly in from off-centre and
settle.

    python scripts/make_logo_svg.py

Animation is plain CSS inside the SVG, because GitHub renders a README image
in an <img> tag, where scripts never run but declarative animation does.
"""

import random
from pathlib import Path

WIDTH, HEIGHT = 560, 120
WORD = "hgt-site"
SUBTITLE = "horizontal gene transfer site caller"
LOOP = 9.0
# Where each letter sits in "hgt-site" at 46px bold in the logo's font,
# measured in a browser rather than guessed, since the font is
# proportional and a fixed step leaves visible gaps.
OFFSETS = [0, 28.63, 57.69, 73.78, 95.30, 121.06, 133.52, 149.92]

# The strand, exactly as the static logo draws it.
HOST_SEGMENTS = [(16, 70), (136, 190)]
INSERT_SEGMENT = (80, 126)
JUNCTIONS = (75, 131)
ROWS = (48, 72)

THEMES = {
    "light": {"host": "#a8a7a1", "insert": "#2a78d6", "word": "#0b0b0b", "sub": "#52514e"},
    "dark": {"host": "#6b6a64", "insert": "#3987e5", "word": "#f0efec", "sub": "#c3c2b7"},
}


def build(theme_name):
    t = THEMES[theme_name]
    rng = random.Random(7)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" role="img" '
        f'aria-label="hgt-site, a horizontal gene transfer site caller">',
        "<style>",
        "g{transform-box:fill-box;transform-origin:center}",
        # Host segments slide in from the outside.
        "@keyframes hostL{0%{opacity:0;transform:translateX(-90px)}"
        "14%,94%{opacity:1;transform:translateX(0)}100%{opacity:0}}",
        "@keyframes hostR{0%{opacity:0;transform:translateX(90px)}"
        "14%,94%{opacity:1;transform:translateX(0)}100%{opacity:0}}",
        # The inserted piece drops in once there is a gap for it.
        "@keyframes drop{0%,13%{opacity:0;transform:translateY(-34px) scaleX(.4)}"
        "26%,94%{opacity:1;transform:translateY(0) scaleX(1)}100%{opacity:0}}",
        # Junction markers appear last, as they do in the analysis.
        "@keyframes mark{0%,27%{opacity:0}36%,94%{opacity:1}100%{opacity:0}}",
        "@keyframes sub{0%,48%{opacity:0}58%,94%{opacity:1}100%{opacity:0}}",
        f".hostL{{animation:hostL {LOOP}s cubic-bezier(.2,.8,.3,1) infinite}}",
        f".hostR{{animation:hostR {LOOP}s cubic-bezier(.2,.8,.3,1) infinite}}",
        f".drop{{animation:drop {LOOP}s cubic-bezier(.2,.8,.3,1) infinite}}",
        f".mark{{animation:mark {LOOP}s ease-out infinite}}",
        f".sub{{animation:sub {LOOP}s ease-out infinite}}",
    ]

    # Letters fly in and settle, staggered, after the strand is together.
    for i in range(len(WORD)):
        delay = 2.9 + 0.11 * i
        dx, dy = rng.uniform(-170, 170), rng.uniform(-80, 80)
        rot = rng.uniform(-140, 140)
        a = delay / LOOP * 100
        b = (delay + 1.0) / LOOP * 100
        out.append(
            f"@keyframes f{i}{{0%,{a:.2f}%{{opacity:0;"
            f"transform:translate({dx:.0f}px,{dy:.0f}px) rotate({rot:.0f}deg) scale(.5)}}"
            f"{b:.2f}%,94%{{opacity:1;transform:translate(0,0) rotate(0) scale(1)}}"
            f"100%{{opacity:0}}}}"
        )
        out.append(f".f{i}{{animation:f{i} {LOOP}s cubic-bezier(.17,.84,.36,1) infinite}}")
    out.append("</style>")

    out.append('<g stroke-linecap="round">')
    for (x1, x2), cls in zip(HOST_SEGMENTS, ("hostL", "hostR")):
        out.append(f'<g class="{cls}">')
        for y in ROWS:
            out.append(f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" stroke="{t["host"]}" stroke-width="6"/>')
        out.append("</g>")
    out.append('<g class="drop">')
    for y in ROWS:
        out.append(
            f'<line x1="{INSERT_SEGMENT[0]}" y1="{y}" x2="{INSERT_SEGMENT[1]}" y2="{y}" '
            f'stroke="{t["insert"]}" stroke-width="6"/>'
        )
    out.append("</g>")
    out.append('<g class="mark">')
    for x in JUNCTIONS:
        out.append(
            f'<line x1="{x}" y1="34" x2="{x}" y2="86" stroke="{t["insert"]}" '
            f'stroke-width="2.5" stroke-dasharray="3 5"/>'
        )
    out.append("</g>")
    out.append("</g>")

    for i, ch in enumerate(WORD):
        out.append(
            f'<g class="f{i}"><text x="{214 + OFFSETS[i]:.2f}" y="68" '
            f'font-family="-apple-system, \'Segoe UI\', Helvetica, Arial, sans-serif" '
            f'font-size="46" font-weight="700" fill="{t["word"]}">{ch}</text></g>'
        )
    out.append(
        f'<text class="sub" x="216" y="96" '
        f'font-family="-apple-system, \'Segoe UI\', Helvetica, Arial, sans-serif" '
        f'font-size="15" fill="{t["sub"]}">{SUBTITLE}</text>'
    )
    out.append("</svg>")
    return "\n".join(out)


if __name__ == "__main__":
    images = Path(__file__).resolve().parent.parent / "docs" / "images"
    for name in THEMES:
        path = images / f"logo-animated-{name}.svg"
        path.write_text(build(name))
        print(f"wrote {path} ({path.stat().st_size:,} bytes)")
