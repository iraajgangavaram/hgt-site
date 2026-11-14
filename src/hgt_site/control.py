"""Mark sites that also appear in a matched control sample.

A tumour and the normal tissue from the same person share everything except
what the tumour did, so a position that produces foreign reads in both of them
is not a somatic integration. It is the protocol, the reference, or something
that was in that person's germline to begin with. That is the strongest
artefact control available for human data, and it needs no model of chimera
rates: the same person's own tissue says what the background looks like.

The control is read from the table hgt-site already writes, so the normal is
called once and the result reused for every tumour compared against it.

What counts as present in the control took measuring rather than reasoning.
The obvious rule is that any candidate in the control counts, however weak,
since an artefact seen ten times in the tumour and twice in the normal is
still an artefact and two reads never pass a filter. That rule is useless. A
control sample produces thousands of weak candidates from chance chimeras,
roughly one per kilobase, so nearly every position in the genome ends up
within reach of one and the control rejects almost everything, real sites
included.

The control has to show the artefact convincingly enough that it is
reproducible rather than noise. A candidate needs five supporting reads to
count, and has to sit within 300 bases of the site. On a tumour and normal
pair sharing a six-copy region, which is the artefact a control exists to
catch, that removed all eight false sites and cost one real one, taking
precision from 0.692 to 1.000.

This does not rescue a host whose artefacts are everywhere. On a bacterial
host, where false sites come from homology with relatives in the sample, the
control carries candidates across the whole genome and no threshold separates
them from the real sites.
"""

import csv
import gzip
from bisect import bisect_left

# How far apart two calls of the same position can sit. A junction estimated
# from read pairs moves by a fragment length between separate runs, so this is
# wider than the clustering window, and no wider: a control sample has about
# one candidate per kilobase, so a window of that size reaches everything.
WINDOW = 300
# Supporting reads the control needs before a site counts as present there.
# Set from the tumour and normal benchmark rather than from first principles;
# see the module docstring.
MIN_PAIRS = 5


def load_control(path, min_pairs=MIN_PAIRS):
    """Positions of the control's candidates, by contig, sorted.

    Every candidate is read, not only the ones that passed. A position held
    back by LOW_SUPPORT in the normal is still evidence about the tumour, so
    long as it has enough reads behind it to be reproducible.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    positions = {}
    with opener(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "chrom" not in reader.fieldnames:
            raise ValueError(
                f"{path} does not look like an hgt-site table: it has no chrom column."
            )
        for row in reader:
            try:
                pairs = int(row.get("n_pairs") or 0)
            except ValueError:
                pairs = 0
            if pairs < min_pairs:
                continue
            start, end = int(row["start"]), int(row["end"])
            positions.setdefault(row["chrom"], []).append((start + end) // 2)
    for chrom in positions:
        positions[chrom].sort()
    return positions


def in_control(positions, chrom, start, end, window=WINDOW):
    """Whether the control has a candidate within `window` of this site."""
    found = positions.get(chrom)
    if not found:
        return False
    # The nearest control position on either side of the site is enough: the
    # list is sorted, so anything closer would be one of these two.
    index = bisect_left(found, start)
    for candidate in found[max(0, index - 1) : index + 1]:
        if candidate >= start - window and candidate <= end + window:
            return True
    return False


def flag_in_control(sites, positions, window=WINDOW):
    """Flag every site the control also has. Returns the count."""
    flagged = 0
    for site in sites:
        if in_control(positions, site.chrom, site.start, site.end, window):
            site.flags.append("IN_CONTROL")
            flagged += 1
    return flagged
