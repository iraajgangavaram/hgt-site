"""Build the animated wordmark shown at the top of the README.

Fragments scatter in from off-centre and settle into the name, which is the
same thing the tool does: loose pieces of sequence resolved into one position.
A few bases drift past behind them.

    python scripts/make_wordmark_svg.py

Animation is plain CSS inside the SVG, because GitHub renders a README image
in an <img> tag, where scripts never run but declarative animation does.
"""

import random
from pathlib import Path

WORD = "hgt-site"
WIDTH, HEIGHT = 760, 150
LETTER_SIZE = 66
BASE_Y = 96
LOOP = 9.0
SETTLE = 1.15

THEMES = {
    "light": {"text": "#23221f", "accent": "#1b7f5a", "drift": "#c9c5bd", "rule": "#dcdad6"},
    "dark": {"text": "#ececea", "accent": "#4fd1a0", "drift": "#3a3f47", "rule": "#2b2f36"},
}


def build(theme_name):
    t = THEMES[theme_name]
    rng = random.Random(11)
    # Monospace advance keeps the letters evenly spaced without measuring.
    advance = LETTER_SIZE * 0.60
    start_x = (WIDTH - advance * len(WORD)) / 2

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" '
        f'font-family="ui-monospace,SFMono-Regular,Menlo,Consolas,monospace" '
        f'role="img" aria-label="hgt-site">',
        "<style>",
        "g{transform-box:fill-box;transform-origin:center}",
    ]

    # Bases drifting across behind the name.
    for i in range(16):
        dur = rng.uniform(5.5, 8.5)
        out.append(
            f"@keyframes d{i}{{0%{{opacity:0;transform:translate(0,0)}}"
            f"18%{{opacity:.55}}82%{{opacity:.55}}"
            f"100%{{opacity:0;transform:translate({rng.uniform(-70, 70):.0f}px,"
            f"{rng.uniform(-26, 26):.0f}px)}}}}"
        )
        out.append(f".d{i}{{animation:d{i} {dur:.2f}s ease-in-out infinite;animation-delay:-{rng.uniform(0, dur):.2f}s}}")

    # Each letter flies in from somewhere else and settles.
    for i in range(len(WORD)):
        delay = 0.13 * i
        dx, dy = rng.uniform(-230, 230), rng.uniform(-120, 120)
        rot = rng.uniform(-150, 150)
        a = delay / LOOP * 100
        b = (delay + SETTLE) / LOOP * 100
        out.append(
            f"@keyframes f{i}{{"
            f"0%,{a:.2f}%{{opacity:0;transform:translate({dx:.0f}px,{dy:.0f}px) rotate({rot:.0f}deg) scale(.55)}}"
            f"{b:.2f}%,92%{{opacity:1;transform:translate(0,0) rotate(0) scale(1)}}"
            f"100%{{opacity:0;transform:translate(0,-10px) scale(.97)}}}}"
        )
        out.append(f".f{i}{{animation:f{i} {LOOP}s cubic-bezier(.17,.84,.36,1) infinite}}")

    # The rule underneath draws itself once the name is together.
    out.append(
        f"@keyframes rule{{0%,16%{{width:0;opacity:0}}34%{{opacity:1}}"
        f"40%,92%{{width:{advance * len(WORD):.0f}px;opacity:1}}100%{{opacity:0}}}}"
    )
    out.append(f".rule{{animation:rule {LOOP}s ease-out infinite}}")
    out.append("</style>")

    for i in range(16):
        x, y = rng.uniform(40, WIDTH - 40), rng.uniform(26, HEIGHT - 20)
        out.append(
            f'<text class="d{i}" x="{x:.0f}" y="{y:.0f}" font-size="{rng.uniform(11, 17):.0f}" '
            f'fill="{t["drift"]}">{rng.choice("ACGT")}</text>'
        )

    for i, ch in enumerate(WORD):
        colour = t["accent"] if ch == "-" else t["text"]
        x = start_x + advance * i
        out.append(
            f'<g class="f{i}"><text x="{x:.1f}" y="{BASE_Y}" font-size="{LETTER_SIZE}" '
            f'font-weight="600" fill="{colour}">{ch}</text></g>'
        )

    out.append(
        f'<rect class="rule" x="{start_x:.0f}" y="{BASE_Y + 17}" height="3" rx="1.5" '
        f'width="0" fill="{t["accent"]}"/>'
    )
    out.append("</svg>")
    return "\n".join(out)


if __name__ == "__main__":
    images = Path(__file__).resolve().parent.parent / "docs" / "images"
    images.mkdir(parents=True, exist_ok=True)
    for name in THEMES:
        path = images / f"wordmark-{name}.svg"
        path.write_text(build(name))
        print(f"wrote {path} ({path.stat().st_size:,} bytes)")
