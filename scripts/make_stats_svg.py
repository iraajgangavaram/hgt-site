"""Build the animated figures panel for the top of the README.

Three numbers that say what the tool is for: how many people have installed it,
how many candidates it has sifted without a false call, and how long it takes.
Each one counts up when the page loads.

    python scripts/make_stats_svg.py

The install count is read from docs/stats.json, which the update-stats workflow
refreshes nightly from the Anaconda API. The other two come from the benchmark
and change only when the benchmark is rerun, so they are kept here.

Counting up is done by drawing every intermediate value as its own <text> and
showing them in sequence, because GitHub renders a README image in an <img>
tag, where scripts never run. CSS keyframes do run, and that is all this needs.
"""

import json
from pathlib import Path

WIDTH, HEIGHT = 900, 116
LOOP = 9.0
# How long the numbers take to land, as a fraction of the loop.
RAMP = 0.34
STEPS = 18

THEMES = {
    "light": {
        "value": "#141312", "label": "#6d6a65", "rule": "#e2dfda",
        "accent": "#157954", "unit": "#8d8a84",
    },
    "dark": {
        "value": "#f0efec", "label": "#9b968e", "rule": "#2b2f36",
        "accent": "#4fd1a0", "unit": "#7e8a97",
    },
}


def ramp(target, steps=STEPS):
    """Intermediate values for a count-up, decelerating towards the target.

    Ease-out rather than linear, so the number slows as it lands instead of
    stopping dead, and no value is repeated.
    """
    out = []
    for i in range(1, steps + 1):
        t = i / steps
        eased = 1 - (1 - t) ** 3
        value = int(round(target * eased))
        if not out or value != out[-1]:
            out.append(value)
    if out[-1] != target:
        out.append(target)
    return out


def counter(x, target, fmt, cls, theme, size=38):
    """A number that counts up to `target` and stays there."""
    values = ramp(target)
    n = len(values)
    out, css = [], []
    for i, value in enumerate(values):
        last = i == n - 1
        # Each frame owns an equal slice of the ramp; the final value holds.
        a = RAMP * i / n * 100
        b = RAMP * (i + 1) / n * 100
        name = f"{cls}{i}"
        if last:
            css.append(f"@keyframes {name}{{0%,{a:.2f}%{{opacity:0}}{min(b, 100):.2f}%,100%{{opacity:1}}}}")
        else:
            css.append(
                f"@keyframes {name}{{0%,{a:.2f}%{{opacity:0}}"
                f"{a + 0.01:.2f}%,{b:.2f}%{{opacity:1}}{b + 0.01:.2f}%,100%{{opacity:0}}}}"
            )
        css.append(f".{name}{{animation:{name} {LOOP}s linear infinite}}")
        out.append(
            f'<text class="{name}" x="{x}" y="58" font-size="{size}" font-weight="700" '
            f'fill="{theme["value"]}" text-anchor="middle" '
            f'font-variant-numeric="tabular-nums">{fmt(value)}</text>'
        )
    return out, css


def build(theme_name, installs):
    t = THEMES[theme_name]
    columns = [
        (installs, lambda v: f"{v:,}", "conda installs", "and counting"),
        (123926, lambda v: f"{v:,}", "candidates sifted", "across 18 benchmark runs"),
        (0, lambda v: str(v), "false positives", "in every one of them"),
    ]
    step = WIDTH / len(columns)
    body, css = [], []

    for i, (target, fmt, label, sub) in enumerate(columns):
        x = step * (i + 0.5)
        texts, rules = counter(x, target, fmt, f"c{i}", t)
        body += texts
        css += rules
        colour = t["accent"] if label == "false positives" else t["label"]
        body.append(
            f'<text x="{x:.0f}" y="84" font-size="14" font-weight="600" fill="{colour}" '
            f'text-anchor="middle">{label}</text>'
        )
        body.append(
            f'<text x="{x:.0f}" y="104" font-size="12" fill="{t["unit"]}" '
            f'text-anchor="middle">{sub}</text>'
        )
        if i:
            body.append(
                f'<line x1="{step * i:.0f}" y1="26" x2="{step * i:.0f}" y2="100" '
                f'stroke="{t["rule"]}" stroke-width="1"/>'
            )

    head = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" '
        f'viewBox="0 0 {WIDTH} {HEIGHT}" '
        f"font-family=\"-apple-system,'Segoe UI',Helvetica,Arial,sans-serif\" role=\"img\" "
        f'aria-label="{installs:,} conda installs. 123,926 candidate sites sifted across 18 '
        f'benchmark runs, with zero false positives in every one of them.">'
    )
    return "\n".join([head, "<style>", *css, "</style>", *body, "</svg>"])


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    stats = json.loads((root / "docs" / "stats.json").read_text())
    installs = int(stats["conda_downloads"])
    for name in THEMES:
        path = root / "docs" / "images" / f"stats-{name}.svg"
        path.write_text(build(name, installs))
        print(f"wrote {path} ({path.stat().st_size:,} bytes)")
