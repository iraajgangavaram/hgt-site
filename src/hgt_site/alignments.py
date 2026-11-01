"""Read host alignments and pull out anchored and fully unmapped read pairs.

BAM and CRAM files are read with pysam. Plain SAM text can also be read
without pysam, which keeps the test suite light.
"""

from collections import Counter
from dataclasses import dataclass

from .kraken import normalise_read_name

PAIRED = 0x1
UNMAPPED = 0x4
MATE_UNMAPPED = 0x8
REVERSE = 0x10
READ1 = 0x40
SECONDARY = 0x100
QCFAIL = 0x200
DUPLICATE = 0x400
SUPPLEMENTARY = 0x800
SKIP = SECONDARY | SUPPLEMENTARY | QCFAIL | DUPLICATE


@dataclass
class Record:
    name: str
    flag: int
    chrom: str
    start: int  # 0-based
    end: int  # 0-based, exclusive
    mapq: int
    seq: str


@dataclass
class Anchor:
    """A host-mapped read whose mate did not map to the host."""

    name: str
    chrom: str
    start: int
    end: int
    reverse: bool
    mapq: int

    @property
    def junction_side(self):
        # A forward anchor points downstream, so the junction is past its end.
        # A reverse anchor points upstream, so the junction is before its start.
        return self.start if self.reverse else self.end


@dataclass
class ReadPairs:
    anchors: list
    unmapped_mate_seqs: dict
    both_unmapped: set
    contig_lengths: dict


def _cigar_ref_length(cigar):
    length, number = 0, ""
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            if char in "MDN=X":
                length += int(number)
            number = ""
    return length


def _wanted(flag):
    """True for primary records of paired reads where at least one mate is unmapped.

    These are the only records hgt-site uses. In a typical host alignment they
    are well under 1% of all records, so checking the flag before building
    anything else saves most of the run time.
    """
    return bool(flag & PAIRED) and not flag & SKIP and bool(flag & (UNMAPPED | MATE_UNMAPPED))


def _iter_sam_text(path):
    lengths = {}
    records = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:])
                    lengths[tags["SN"]] = int(tags["LN"])
                continue
            cols = line.split("\t", 10)
            flag = int(cols[1])
            if not _wanted(flag):
                continue
            pos = int(cols[3]) - 1
            cigar = cols[5]
            end = pos + (_cigar_ref_length(cigar) if cigar != "*" else 0)
            records.append(Record(cols[0], flag, cols[2], pos, end, int(cols[4]), cols[9]))
    return lengths, records


def _iter_pysam(path, reference=None, threads=1):
    import pysam

    mode = "rc" if str(path).endswith(".cram") else "rb"
    bam = pysam.AlignmentFile(str(path), mode, reference_filename=reference, threads=threads)
    lengths = dict(zip(bam.references, bam.lengths))

    def generate():
        for read in bam.fetch(until_eof=True):
            if not _wanted(read.flag):
                continue
            yield Record(
                read.query_name,
                read.flag,
                read.reference_name or "*",
                read.reference_start if read.reference_start is not None else -1,
                read.reference_end if read.reference_end is not None else -1,
                read.mapping_quality,
                read.query_sequence or "",
            )

    return lengths, generate()


def open_alignments(path, reference=None, threads=1):
    path = str(path)
    if path.endswith(".sam"):
        return _iter_sam_text(path)
    return _iter_pysam(path, reference, threads)


def collect_pairs(path, min_mapq=20, reference=None, threads=1):
    """Scan the alignments once and sort read pairs into the groups we need."""
    lengths, records = open_alignments(path, reference, threads)
    anchors = []
    mate_seqs = {}
    both_unmapped = set()
    for rec in records:
        if not rec.flag & PAIRED or rec.flag & SKIP:
            continue
        name = normalise_read_name(rec.name)
        read_unmapped = bool(rec.flag & UNMAPPED)
        mate_unmapped = bool(rec.flag & MATE_UNMAPPED)
        if read_unmapped and mate_unmapped:
            both_unmapped.add(name)
        elif not read_unmapped and mate_unmapped:
            if rec.mapq >= min_mapq:
                anchors.append(
                    Anchor(name, rec.chrom, rec.start, rec.end, bool(rec.flag & REVERSE), rec.mapq)
                )
        elif read_unmapped and not mate_unmapped:
            mate_seqs[name] = rec.seq
    return ReadPairs(anchors, mate_seqs, both_unmapped, lengths)


def count_by_taxon(names, calls, taxon_filter):
    counts = Counter()
    for name in names:
        taxid = calls.get(name)
        if taxid is not None and taxon_filter.is_foreign(taxid):
            counts[taxid] += 1
    return counts
