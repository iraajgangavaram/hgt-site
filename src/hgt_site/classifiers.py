"""Read per-read taxonomic calls from whichever classifier the user ran.

hgt-site needs one thing from a classifier: which taxon each unmapped read was
assigned to. Requiring that to come from Kraken2 specifically excluded everyone
who classified their reads another way, which was a limit of the tool rather
than of the method.

Supported, and detected from the file itself:

- Kraken2 and KrakenUniq per-read output, and Kaiju, which share a layout:
  C or U, read name, taxid.
- Centrifuge, which writes a header line starting readID and puts the taxid in
  the third column.
- DIAMOND with --outfmt 102: read name, taxid, e-value.
- BLAST tabular output (-outfmt 6) ending in staxids.
- Any two-column file of read name and taxid, which is the fallback for a
  classifier nobody has thought of.

For formats that report several hits per read, the first line for a read wins,
because BLAST and DIAMOND both write the best hit first.
"""

import gzip
import logging
import re
from pathlib import Path

log = logging.getLogger("hgt_site")

_TAXID_IN_NAME = re.compile(r"\(taxid (\d+)\)")

FORMATS = ("auto", "kraken2", "kaiju", "centrifuge", "diamond", "blast", "tsv")


def normalise_read_name(name):
    """Drop the /1 or /2 mate suffix so both mates share one name."""
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def parse_taxid(field):
    """Read a taxid from a field, which may be a bare number or a name with one.

    Several taxids separated by ; or , means the classifier could not decide
    between them, and the first is taken.
    """
    field = field.strip()
    if ";" in field or "," in field:
        field = re.split("[;,]", field)[0].strip()
    if field.isdigit():
        return int(field)
    match = _TAXID_IN_NAME.search(field)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not read a taxid from the field {field!r}")


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _first_data_line(path):
    with _open(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\n")
    return ""


def detect_format(path):
    """Work out which classifier wrote a file, from its first line of data."""
    line = _first_data_line(path)
    if not line:
        return "tsv"
    cols = line.split("\t")
    if cols[0] == "readID" or (len(cols) >= 3 and cols[2] == "taxID"):
        return "centrifuge"
    if cols[0] in ("C", "U") and len(cols) >= 3:
        return "kraken2"
    if len(cols) == 2:
        return "tsv"
    if len(cols) >= 12:
        return "blast"
    if len(cols) == 3:
        return "diamond"
    raise ValueError(
        f"Could not tell which classifier wrote {path}. Its first line has {len(cols)} "
        "tab-separated columns, which matches none of the formats hgt-site knows. "
        "Pass --classifier to say which it is, or give a two-column file of read name "
        "and taxid."
    )


def _rows(path, fmt):
    """Yield (read name, taxid field) pairs for one file in one format."""
    with _open(path) as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            try:
                if fmt in ("kraken2", "kaiju"):
                    if cols[0] != "C":
                        continue
                    yield cols[1], cols[2]
                elif fmt == "centrifuge":
                    if cols[0] == "readID":
                        continue
                    yield cols[0], cols[2]
                elif fmt == "diamond":
                    yield cols[0], cols[1]
                elif fmt == "blast":
                    yield cols[0], cols[-1]
                else:
                    yield cols[0], cols[1]
            except IndexError:
                raise ValueError(
                    f"{path} line {number} has {len(cols)} columns, too few for {fmt} output."
                ) from None


def read_calls(path, fmt="auto"):
    """Return {read name: taxid} for one file of per-read classifications."""
    if fmt == "auto":
        fmt = detect_format(path)
    calls = {}
    unreadable = 0
    for name, field in _rows(path, fmt):
        name = normalise_read_name(name)
        if name in calls:
            # The first line for a read wins: BLAST and DIAMOND write the best
            # hit first, and a classifier that reports a read once is unaffected.
            continue
        try:
            taxid = parse_taxid(field)
        except ValueError:
            unreadable += 1
            continue
        calls[name] = taxid
    if unreadable:
        log.warning("Skipped %d lines of %s whose taxid could not be read.", unreadable, path)
    if not calls:
        log.warning("No classified reads were found in %s, read as %s output.", path, fmt)
    return calls, fmt


def load_classifications(paths, fmt="auto"):
    """Read one or more per-read classification files into one {read: taxid} map."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    calls = {}
    seen = set()
    for path in paths:
        part, used = read_calls(Path(path), fmt)
        seen.add(used)
        calls.update(part)
    if calls:
        log.info(
            "Read %d classified reads from %d file(s) as %s output.",
            len(calls), len(paths), " and ".join(sorted(seen)),
        )
    return calls
