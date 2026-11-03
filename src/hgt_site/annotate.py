"""Say which gene a candidate site falls in, or which gene is nearest.

A position on its own is not much use to a biologist. The first question
anyone asks about a candidate insertion is which gene it landed in, and the
second is what the nearest gene is when it landed between genes. This module
answers both from a GTF or GFF3 annotation, which is the file people already
have beside their reference genome.

Only gene records are read. Transcript and exon records would let a site be
placed in an intron or an exon, but they multiply the file size by ten and
the gene is what gets reported in a paper.
"""

import bisect
import gzip
import logging
import re

log = logging.getLogger("hgt_site")

# How far from a gene a site can be and still have that gene reported as the
# nearest one. Beyond this the site is called intergenic and no gene is given,
# because a gene half a megabase away says nothing about the insertion.
MAX_DISTANCE = 100_000

_GTF_NAME = re.compile(r'gene_name "([^"]+)"')
_GTF_ID = re.compile(r'gene_id "([^"]+)"')
_GFF_NAME = re.compile(r"(?:^|;)(?:Name|gene_name)=([^;]+)")
_GFF_ID = re.compile(r"(?:^|;)(?:ID|gene_id)=([^;]+)")


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def _names(attributes):
    """Pull the gene name and identifier out of a GTF or GFF3 attribute field."""
    name = _GTF_NAME.search(attributes) or _GFF_NAME.search(attributes)
    identifier = _GTF_ID.search(attributes) or _GFF_ID.search(attributes)
    name = name.group(1) if name else ""
    identifier = identifier.group(1) if identifier else ""
    if identifier.startswith("gene:"):
        # Ensembl GFF3 writes ID=gene:ENSG00000141510.
        identifier = identifier[len("gene:"):]
    return name or identifier, identifier or name


def load_genes(path):
    """Read gene records from a GTF or GFF3 file, gzipped or not.

    Returns {chrom: (starts, records)} where starts is sorted for bisect and
    records holds (start, end, name, gene_id) with 0-based half-open
    coordinates, matching the rest of hgt-site.
    """
    by_chrom = {}
    kept = 0
    with _open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            try:
                start, end = int(cols[3]) - 1, int(cols[4])
            except ValueError:
                continue
            name, identifier = _names(cols[8])
            if not name:
                continue
            by_chrom.setdefault(cols[0], []).append((start, end, name, identifier))
            kept += 1
    if not kept:
        raise ValueError(
            f"No gene records were found in {path}. hgt-site reads the lines whose third "
            "column is 'gene', which GTF and GFF3 files from Ensembl, GENCODE and NCBI all have."
        )
    indexed = {}
    for chrom, records in by_chrom.items():
        records.sort()
        indexed[chrom] = ([r[0] for r in records], records)
    log.info("Read %d genes on %d sequences from the annotation.", kept, len(indexed))
    return indexed


def _nearest(genes, chrom, position, max_distance):
    """Gene containing position, else the closest one within max_distance.

    Returns (name, gene_id, distance) with distance 0 when the position is
    inside the gene, or None when nothing is close enough.
    """
    entry = genes.get(chrom)
    if entry is None:
        return None
    starts, records = entry
    # Genes starting after the position cannot contain it, so walk back from
    # the first one that does. Overlapping genes mean the containing gene is
    # not always the last one to start, hence the window rather than one step.
    index = bisect.bisect_right(starts, position)
    best = None
    for record in records[max(0, index - 50):index]:
        start, end, name, identifier = record
        if start <= position < end:
            return name, identifier, 0
        distance = position - end + 1
        if distance >= 0 and (best is None or distance < best[2]):
            best = (name, identifier, distance)
    for record in records[index:index + 50]:
        start, end, name, identifier = record
        distance = start - position
        if best is None or distance < best[2]:
            best = (name, identifier, distance)
        break
    if best is None or best[2] > max_distance:
        return None
    return best


def annotate(sites, path, junction, max_distance=MAX_DISTANCE):
    """Add gene, gene_id and gene_distance to every site's stats.

    junction is the function giving a site's 1-based reported position, so the
    gene reported is the gene at the position hgt-site actually calls rather
    than at the edge of the clustered region.
    """
    genes = load_genes(path)
    found = 0
    for site in sites:
        hit = _nearest(genes, site.chrom, junction(site) - 1, max_distance)
        if hit is None:
            site.stats.update(gene="", gene_id="", gene_distance="")
        else:
            name, identifier, distance = hit
            site.stats.update(gene=name, gene_id=identifier, gene_distance=distance)
            found += 1
    return found
