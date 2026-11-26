"""Build a portable, browser-readable report from an hgt-site result table."""

import csv
import html
import json
import os
from collections import Counter
from pathlib import Path
from urllib.parse import quote


REQUIRED_COLUMNS = {"chrom", "start", "end", "n_pairs", "top_taxid", "filter"}
DISPLAY_COLUMNS = (
    "chrom",
    "start",
    "end",
    "gene",
    "top_taxid",
    "n_pairs",
    "qvalue",
    "left_breakpoint",
    "left_split_reads",
    "right_breakpoint",
    "right_split_reads",
    "filter",
    "flags",
)
TABLE_COLUMNS = (
    ("chrom", "Chromosome"),
    ("start", "Start"),
    ("end", "End"),
    ("gene", "Gene"),
    ("top_taxid", "Taxid"),
    ("n_pairs", "Pairs"),
    ("qvalue", "Q-value"),
    ("junctions", "Junctions"),
    ("filter", "Filter"),
    ("flags", "Flags"),
)


def _read_rows(path):
    """Read the call table and reject rows that cannot be shown faithfully."""
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} is not an hgt-site table: no header row.")
        duplicates = sorted({name for name in reader.fieldnames if reader.fieldnames.count(name) > 1})
        if duplicates:
            raise ValueError(f"{path} has duplicate column(s): {', '.join(duplicates)}.")
        missing = sorted(REQUIRED_COLUMNS - set(reader.fieldnames))
        if missing:
            raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}.")

        for row in reader:
            line = reader.line_num
            absent = [name for name in reader.fieldnames if row.get(name) is None]
            if absent or None in row:
                problems = []
                if absent:
                    problems.append(f"missing value(s) for {', '.join(absent)}")
                if None in row:
                    problems.append(f"{len(row[None])} extra field(s)")
                raise ValueError(f"{path}:{line}: malformed row ({'; '.join(problems)}).")
            rows.append({name: row.get(name, "") or "" for name in DISPLAY_COLUMNS})
    return rows


def _summary(rows):
    passed = sum(row["filter"] == "PASS" for row in rows)
    exact = sum(
        any(value and value != "." for value in (row["left_breakpoint"], row["right_breakpoint"]))
        for row in rows
    )
    failure_reasons = Counter()
    taxa = Counter()
    for row in rows:
        taxid = row["top_taxid"]
        if taxid and taxid != ".":
            taxa[taxid] += 1
        if row["filter"] != "PASS":
            failure_reasons.update(
                flag for flag in row["filter"].split(";") if flag and flag != "."
            )
    return passed, exact, failure_reasons, taxa


def _table_records(rows):
    records = []
    for row in rows:
        left = row["left_breakpoint"]
        right = row["right_breakpoint"]
        junctions = []
        if left:
            junctions.append(f"L: {left} ({row['left_split_reads']} reads)")
        if right:
            junctions.append(f"R: {right} ({row['right_split_reads']} reads)")
        record = dict(row)
        record["junctions"] = "; ".join(junctions) or "Estimated from read pairs"
        records.append(record)
    return records


def render_report(rows, source_name, source_href):
    """Return a self-contained HTML report for parsed hgt-site rows."""
    total = len(rows)
    passed, exact, reasons, taxa = _summary(rows)
    records = _table_records(rows)
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":")).replace("<", "\\u003c")
    title = f"HGT-Site report: {source_name}"

    def summary_table(counter, empty):
        if not counter:
            return f"<p class=\"muted\">{empty}</p>"
        body = "".join(
            f"<tr><th scope=\"row\">{html.escape(str(name))}</th>"
            f"<td>{count:,}</td></tr>"
            for name, count in counter.most_common(10)
        )
        return f"<table class=\"summary-table\"><tbody>{body}</tbody></table>"

    headings = "".join(f"<th scope=\"col\">{html.escape(label)}</th>" for _, label in TABLE_COLUMNS)
    source = html.escape(source_name)
    source_link = html.escape(source_href, quote=True)
    escaped_title = html.escape(title)
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escaped_title}</title>
<style>
:root {{ color-scheme: light; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #17212b; background: #f4f7fa; }}
body {{ margin: 0; }}
main {{ max-width: 1200px; margin: 0 auto; padding: 2rem 1.25rem 4rem; }}
h1 {{ margin: 0 0 .4rem; font-size: clamp(1.7rem, 4vw, 2.5rem); }}
h2 {{ margin: 2rem 0 .7rem; font-size: 1.25rem; }}
p {{ line-height: 1.55; }}
.muted {{ color: #5b6875; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(145px, 1fr)); gap: .8rem; margin: 1.4rem 0; }}
.card {{ background: white; border: 1px solid #d9e1e8; border-radius: .6rem; padding: 1rem; }}
.card strong {{ display: block; font-size: 1.8rem; }}
.card span {{ color: #52606d; font-size: .9rem; }}
.summaries {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 1.5rem; }}
.summary-table {{ border-collapse: collapse; width: 100%; background: white; }}
.summary-table th, .summary-table td {{ border-bottom: 1px solid #e1e7ec; padding: .55rem .65rem; text-align: left; }}
.summary-table td {{ text-align: right; font-variant-numeric: tabular-nums; }}
.controls {{ display: flex; flex-wrap: wrap; gap: .7rem; align-items: center; margin: .8rem 0; }}
input, select, button {{ font: inherit; border: 1px solid #aebbc7; border-radius: .35rem; padding: .5rem .65rem; background: white; }}
input {{ flex: 1 1 18rem; }}
button {{ cursor: pointer; }}
.table-wrap {{ overflow-x: auto; background: white; border: 1px solid #d9e1e8; border-radius: .5rem; }}
.calls {{ border-collapse: collapse; width: 100%; min-width: 900px; }}
.calls th, .calls td {{ border-bottom: 1px solid #e4e9ee; padding: .55rem .65rem; text-align: left; vertical-align: top; }}
.calls thead {{ background: #eaf1f7; position: sticky; top: 0; }}
.calls td:nth-child(2), .calls td:nth-child(3), .calls td:nth-child(6), .calls td:nth-child(7) {{ font-variant-numeric: tabular-nums; }}
.status-pass {{ color: #176b3a; font-weight: 650; }}
.status-fail {{ color: #8c3a23; }}
.pager {{ display: flex; gap: .6rem; align-items: center; justify-content: flex-end; margin: .7rem 0; }}
@media print {{ :root {{ background: white; }} main {{ max-width: none; padding: 0; }} .controls, .pager {{ display: none; }} .table-wrap {{ overflow: visible; }} .calls {{ min-width: 0; }} .calls thead {{ position: static; }} }}
</style>
</head>
<body>
<main>
<h1>{escaped_title}</h1>
<p class="muted">Summary of <a href="{source_link}">{source}</a>. This report does not rerun analysis or change the caller's filters.</p>
<section class="cards" aria-label="Call summary">
  <div class="card"><strong>{total:,}</strong><span>candidate sites</span></div>
  <div class="card"><strong>{passed:,}</strong><span>passed filters</span></div>
  <div class="card"><strong>{total - passed:,}</strong><span>did not pass</span></div>
  <div class="card"><strong>{exact:,}</strong><span>with split-read junctions</span></div>
</section>
<section class="summaries">
  <div><h2>Filter reasons</h2>{summary_table(reasons, "No candidates failed a filter.")}</div>
  <div><h2>Taxids in candidate calls</h2>{summary_table(taxa, "No taxids are recorded in this table.")}</div>
</section>
<section>
<h2>Candidate sites</h2>
<p class="muted">Search by chromosome, gene, taxid, filter, or any visible value. Coordinates are shown as recorded in the TSV.</p>
<div class="controls">
  <label for="search">Search</label><input id="search" type="search" placeholder="e.g. chr19, PASS, or gene name">
  <label for="filter">Show</label><select id="filter"><option value="all">All candidates</option><option value="PASS">Passed only</option><option value="failed">Did not pass</option></select>
  <button type="button" onclick="window.print()">Print report</button>
</div>
<div class="table-wrap"><table class="calls"><thead><tr>{headings}</tr></thead><tbody id="rows"></tbody></table></div>
<div class="pager"><span id="count" aria-live="polite"></span><button id="previous" type="button">Previous</button><button id="next" type="button">Next</button></div>
</section>
<noscript><p>Enable JavaScript in your browser to search and page through the candidate table.</p></noscript>
</main>
<script>
const DATA = {payload};
const COLUMNS = {json.dumps([key for key, _ in TABLE_COLUMNS])};
const PAGE_SIZE = 50;
let currentPage = 0;
function escapeCell(value) {{
  return String(value).replace(/[&<>"']/g, char => ({{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}}[char]));
}}
function visibleRows() {{
  const query = document.getElementById("search").value.toLowerCase().trim();
  const state = document.getElementById("filter").value;
  return DATA.filter(row => {{
    const isPass = row.filter === "PASS";
    return (state === "all" || (state === "PASS" && isPass) || (state === "failed" && !isPass)) &&
      (!query || COLUMNS.some(key => String(row[key] || "").toLowerCase().includes(query)));
  }});
}}
function render(printAll = false) {{
  const rows = visibleRows();
  const pages = Math.max(1, Math.ceil(rows.length / PAGE_SIZE));
  currentPage = Math.min(currentPage, pages - 1);
  const start = currentPage * PAGE_SIZE;
  const body = document.getElementById("rows");
  const pageRows = printAll ? rows : rows.slice(start, start + PAGE_SIZE);
  body.innerHTML = pageRows.map(row => {{
    const values = COLUMNS.map(key => `<td>${{escapeCell(row[key] || "")}}</td>`);
    const statusIndex = COLUMNS.indexOf("filter");
    values[statusIndex] = `<td class="${{row.filter === "PASS" ? "status-pass" : "status-fail"}}">${{escapeCell(row.filter || "")}}</td>`;
    return `<tr>${{values.join("")}}</tr>`;
  }}).join("");
  document.getElementById("count").textContent = rows.length
    ? (printAll ? `Printing all ${{rows.length}} matching candidates` : `Showing ${{start + 1}}–${{Math.min(start + PAGE_SIZE, rows.length)}} of ${{rows.length}}`)
    : "No candidate sites match this search.";
  document.getElementById("previous").disabled = currentPage === 0;
  document.getElementById("next").disabled = currentPage >= pages - 1;
}}
document.getElementById("search").addEventListener("input", () => {{ currentPage = 0; render(); }});
document.getElementById("filter").addEventListener("change", () => {{ currentPage = 0; render(); }});
document.getElementById("previous").addEventListener("click", () => {{ currentPage--; render(); }});
document.getElementById("next").addEventListener("click", () => {{ currentPage++; render(); }});
window.addEventListener("beforeprint", () => render(true));
window.addEventListener("afterprint", () => render());
render();
</script>
</body>
</html>
"""


def write_report(sites_path, output_path):
    """Write a self-contained HTML report and return the number of candidates."""
    rows = _read_rows(sites_path)
    output = Path(output_path)
    source = Path(sites_path)
    source_href = quote(Path(os.path.relpath(source, output.parent)).as_posix(), safe="/")
    document = render_report(rows, source.name, source_href)
    with open(output, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(document)
    return len(rows)
