"""Read Kraken2 per-read output and decide which taxa count as non-host."""

import re
from pathlib import Path

_TAXID_IN_NAME = re.compile(r"\(taxid (\d+)\)")

# Taxids that tell us nothing about where a read came from:
# 0 = unclassified, 1 = root, 131567 = cellular organisms.
UNINFORMATIVE = {0, 1, 131567}


def normalise_read_name(name):
    """Drop the /1 or /2 mate suffix so both mates share one name."""
    if name.endswith("/1") or name.endswith("/2"):
        return name[:-2]
    return name


def parse_taxid(field):
    field = field.strip()
    if field.isdigit():
        return int(field)
    match = _TAXID_IN_NAME.search(field)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not read a taxid from Kraken2 field: {field!r}")


def read_kraken(path):
    """Return {read name: taxid} for every classified read in a Kraken2 output file."""
    calls = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                raise ValueError(f"Not a Kraken2 per-read line: {line!r}")
            if cols[0] != "C":
                continue
            calls[normalise_read_name(cols[1])] = parse_taxid(cols[2])
    return calls


def read_nodes(path):
    """Read an NCBI nodes.dmp file into a {taxid: parent taxid} map."""
    parents = {}
    with open(path) as handle:
        for line in handle:
            cols = [c.strip() for c in line.split("|")]
            parents[int(cols[0])] = int(cols[1])
    return parents


class TaxonFilter:
    """Decides whether a taxid is informative and outside the host lineage.

    Without a taxonomy, only the host taxid itself and the uninformative
    taxids are excluded. With nodes.dmp, anything inside the host clade and
    every ancestor of the host (for example Eukaryota) is excluded too.
    """

    def __init__(self, host_taxid, nodes_dmp=None, exclude=()):
        self.host = int(host_taxid)
        self.extra = {int(t) for t in exclude}
        self.parents = read_nodes(nodes_dmp) if nodes_dmp else None
        self._cache = {}
        self.host_ancestors = self._lineage(self.host) if self.parents else {self.host}

    def _lineage(self, taxid):
        seen = []
        current = taxid
        while current not in seen:
            seen.append(current)
            parent = self.parents.get(current)
            if parent is None or parent == current:
                break
            current = parent
        return set(seen)

    def is_foreign(self, taxid):
        if taxid in self._cache:
            return self._cache[taxid]
        if taxid in UNINFORMATIVE or taxid in self.extra or taxid == self.host:
            result = False
        elif self.parents is None:
            result = True
        else:
            lineage = self._lineage(taxid)
            in_host_clade = self.host in lineage
            is_host_ancestor = taxid in self.host_ancestors
            result = not (in_host_clade or is_host_ancestor)
        self._cache[taxid] = result
        return result


def load_classifications(paths):
    """Read one or more Kraken2 per-read files into a single {read name: taxid} map."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    calls = {}
    for path in paths:
        calls.update(read_kraken(Path(path)))
    return calls
