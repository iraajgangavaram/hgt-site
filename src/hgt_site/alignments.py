"""Read host alignments and pull out anchored and fully unmapped read pairs.

BAM and CRAM files are read with pysam. Plain SAM text can also be read
without pysam, which keeps the test suite light.
"""

from collections import Counter
from dataclasses import dataclass, field

from .composition import Background
from .kraken import normalise_read_name

PAIRED = 0x1
PROPER_PAIR = 0x2
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


# Width of the bins host coverage is counted in. Wide enough that a bin holds
# enough reads for a steady estimate, narrow enough to follow an amplified
# region or a highly expressed gene.
COVERAGE_BIN = 10_000


class Coverage:
    """Mapped reads per bin, for estimating local host coverage.

    Counted inside the readers rather than in collect_pairs, because the
    readers drop the fully mapped pairs before anything else is done with
    them and those pairs are most of the host coverage. Only built when the
    caller asks for it, so the usual path keeps its speed.
    """

    __slots__ = ("bins", "bin_size", "mapped_reads", "fully_mapped_reads")

    def __init__(self, bin_size=None, bins=None, mapped_reads=0, fully_mapped_reads=0):
        self.bins = Counter() if bins is None else bins
        self.bin_size = bin_size
        self.mapped_reads = mapped_reads
        self.fully_mapped_reads = fully_mapped_reads

    def add(self, flag, chrom, pos):
        if flag & UNMAPPED or pos < 0 or chrom == "*":
            return
        self.bins[(chrom, pos // self.bin_size)] += 1
        self.mapped_reads += 1
        if not flag & MATE_UNMAPPED:
            self.fully_mapped_reads += 1

    def rate(self, chrom, start, end, contig_length=None):
        """Mapped reads per base over [start, end), from the bins it touches.

        The whole of every bin the range touches is used, so the estimate is
        over at least one bin's width however narrow the range is. A site is a
        few hundred bases, where read counts are far too noisy to divide by.

        The last bin of a contig stops at the contig end rather than a bin
        boundary, so its width is clipped. Without that, a site near the end
        of a contig would be divided by more bases than exist there, giving a
        rate below the truth and an expected count too low, which is the
        direction that lets a false site through.
        """
        first, last = start // self.bin_size, (max(end - 1, start)) // self.bin_size
        reads = sum(self.bins.get((chrom, b), 0) for b in range(first, last + 1))
        width = (last - first + 1) * self.bin_size
        if contig_length is not None:
            width = min(width, max(contig_length - first * self.bin_size, 1))
        return reads / width, reads

    def mean_rate(self, genome_length):
        return self.mapped_reads / genome_length if genome_length else 0.0


@dataclass
class ReadPairs:
    anchors: list
    unmapped_mate_seqs: dict
    both_unmapped: set
    contig_lengths: dict
    insert_sizes: list = field(default_factory=list)
    read_lengths: list = field(default_factory=list)
    # Sequences of the pairs with both mates unmapped, kept only when the
    # caller asks, since on a whole-genome file there can be a great many.
    both_unmapped_seqs: dict = field(default_factory=dict)
    coverage: "Coverage" = None
    host_composition: "Background" = None


# Number of properly paired reads sampled to estimate the fragment size.
FRAGMENT_SAMPLE = 100_000


class _Sample:
    """Insert sizes and read lengths from the first properly paired reads seen."""

    def __init__(self, limit=FRAGMENT_SAMPLE):
        self.limit = limit
        self.insert_sizes = []
        self.read_lengths = []

    def add(self, flag, tlen, read_length):
        """Add one record; return False once the sample is full."""
        if (flag & PROPER_PAIR) and (flag & READ1) and not (flag & SKIP) and tlen > 0:
            self.insert_sizes.append(tlen)
            self.read_lengths.append(read_length)
        return len(self.insert_sizes) < self.limit


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


def _iter_sam_text(path, sample, coverage=None):
    lengths = {}
    records = []
    sampling = True
    with open(path) as handle:
        for line in handle:
            if line.startswith("@"):
                if line.startswith("@SQ"):
                    tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:])
                    lengths[tags["SN"]] = int(tags["LN"])
                continue
            cols = line.split("\t", 10)
            flag = int(cols[1])
            if sampling:
                sampling = sample.add(flag, int(cols[8]), len(cols[9]))
            if coverage is not None:
                coverage.add(flag, cols[2], int(cols[3]) - 1)
            if not _wanted(flag):
                continue
            pos = int(cols[3]) - 1
            cigar = cols[5]
            end = pos + (_cigar_ref_length(cigar) if cigar != "*" else 0)
            records.append(Record(cols[0], flag, cols[2], pos, end, int(cols[4]), cols[9]))
    return lengths, records


def _iter_pysam(path, sample, reference=None, threads=1, coverage=None):
    import pysam

    mode = "rc" if str(path).endswith(".cram") else "rb"
    bam = pysam.AlignmentFile(str(path), mode, reference_filename=reference, threads=threads)
    lengths = dict(zip(bam.references, bam.lengths))

    def generate():
        sampling = True
        for read in bam.fetch(until_eof=True):
            flag = read.flag
            if sampling:
                sampling = sample.add(flag, read.template_length, read.query_length)
            if coverage is not None:
                coverage.add(
                    flag, read.reference_name or "*",
                    read.reference_start if read.reference_start is not None else -1,
                )
            if not _wanted(flag):
                continue
            yield Record(
                read.query_name,
                flag,
                read.reference_name or "*",
                read.reference_start if read.reference_start is not None else -1,
                read.reference_end if read.reference_end is not None else -1,
                read.mapping_quality,
                read.query_sequence or "",
            )

    return lengths, generate()


def open_alignments(path, sample, reference=None, threads=1, coverage=None):
    path = str(path)
    if path.endswith(".sam"):
        return _iter_sam_text(path, sample, coverage)
    return _iter_pysam(path, sample, reference, threads, coverage)


def collect_pairs(
    path, min_mapq=20, reference=None, threads=1, keep_both_unmapped_seqs=False,
    count_coverage=False, bin_size=COVERAGE_BIN,
):
    """Scan the alignments once and sort read pairs into the groups we need."""
    sample = _Sample()
    coverage = Coverage(bin_size) if count_coverage else None
    host_composition = Background()
    lengths, records = open_alignments(path, sample, reference, threads, coverage)
    anchors = []
    mate_seqs = {}
    both_unmapped = set()
    both_unmapped_seqs = {}
    for rec in records:
        if not rec.flag & PAIRED or rec.flag & SKIP:
            continue
        name = normalise_read_name(rec.name)
        read_unmapped = bool(rec.flag & UNMAPPED)
        mate_unmapped = bool(rec.flag & MATE_UNMAPPED)
        if read_unmapped and mate_unmapped:
            both_unmapped.add(name)
            if keep_both_unmapped_seqs and name not in both_unmapped_seqs:
                both_unmapped_seqs[name] = rec.seq
        elif not read_unmapped and mate_unmapped:
            if rec.mapq >= min_mapq:
                anchors.append(
                    Anchor(name, rec.chrom, rec.start, rec.end, bool(rec.flag & REVERSE), rec.mapq)
                )
                # This read is host sequence by definition, so it gives the
                # host's own base composition for nothing.
                host_composition.add(rec.seq)
        elif read_unmapped and not mate_unmapped:
            mate_seqs[name] = rec.seq
    return ReadPairs(
        anchors, mate_seqs, both_unmapped, lengths, sample.insert_sizes, sample.read_lengths,
        both_unmapped_seqs, coverage, host_composition,
    )


def count_by_taxon(names, calls, taxon_filter):
    counts = Counter()
    for name in names:
        taxid = calls.get(name)
        if taxid is not None and taxon_filter.is_foreign(taxid):
            counts[taxid] += 1
    return counts
