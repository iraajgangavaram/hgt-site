"""Decide which taxa count as non-host, and read per-read classifications.

The reading of classifier output lives in classifiers.py, which handles the
several formats people arrive with. This module keeps the taxonomy side and
re-exports the readers so that existing code and scripts keep working.
"""

from .classifiers import (  # noqa: F401
    detect_format,
    load_classifications,
    normalise_read_name,
    parse_taxid,
    read_calls,
)

# Taxids that tell us nothing about where a read came from:
# 0 = unclassified, 1 = root, 131567 = cellular organisms.
UNINFORMATIVE = {0, 1, 131567}


def read_kraken(path):
    """Kept for callers that used it directly; reads any supported format."""
    calls, _ = read_calls(path)
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
