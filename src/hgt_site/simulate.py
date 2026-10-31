"""Simulate a host sample with known foreign insertions, for testing and benchmarking.

The simulator writes the files hgt-site reads (a SAM file of reads "aligned"
to the host and a Kraken2-style per-read file) plus a BED file of the true
insertion sites. Alignment and classification are modelled directly rather
than run with real tools, so the truth is known exactly. It also adds the two
kinds of noise that make real data hard:

* free contaminant read pairs, where neither mate maps to the host
* chimeric artefact pairs, where a host fragment and a contaminant fragment
  were joined during library preparation at a random host position
"""

import random
from dataclasses import dataclass

HOST_TAXID = 9606
DONOR_TAXID = 562
CONTAMINANT_TAXID = 1280
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def random_seq(rng, length, gc=0.5):
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


@dataclass
class Segment:
    start: int  # in the sample chromosome
    end: int
    source: str  # "host" or a taxid
    src_start: int  # host coordinate for host segments


@dataclass
class SimConfig:
    n_chroms: int = 2
    chrom_length: int = 100_000
    n_insertions: int = 5
    insert_min: int = 1_500
    insert_max: int = 3_000
    coverage: float = 10.0
    read_length: int = 100
    fragment_mean: int = 350
    fragment_sd: int = 30
    contaminant_pairs: int = 2_000
    chimera_pairs: int = 150
    kraken_sensitivity: float = 0.9
    min_mapped_fraction: float = 0.7
    seed: int = 1


def _build_sample_chrom(rng, host, donor, n_insert, cfg, chrom):
    positions = sorted(rng.sample(range(5_000, len(host) - 5_000), n_insert))
    segments, pieces, truth = [], [], []
    cursor = sample_pos = 0
    for pos in positions:
        length = rng.randint(cfg.insert_min, cfg.insert_max)
        dstart = rng.randint(0, len(donor) - length)
        segments.append(Segment(sample_pos, sample_pos + pos - cursor, "host", cursor))
        pieces.append(host[cursor:pos])
        sample_pos += pos - cursor
        segments.append(Segment(sample_pos, sample_pos + length, str(DONOR_TAXID), -1))
        pieces.append(donor[dstart : dstart + length])
        sample_pos += length
        cursor = pos
        truth.append((chrom, pos, length))
    segments.append(Segment(sample_pos, sample_pos + len(host) - cursor, "host", cursor))
    pieces.append(host[cursor:])
    return "".join(pieces), segments, truth


def _place(segments, start, length, cfg):
    """Decide how a read at sample coordinates [start, start+length) aligns to the host.

    Returns (mapped, host_start, cigar, foreign_taxid).
    """
    end = start + length
    host_bases, foreign = [], {}
    for seg in segments:
        lo, hi = max(start, seg.start), min(end, seg.end)
        if lo >= hi:
            continue
        if seg.source == "host":
            host_bases.append((lo, hi, seg))
        else:
            foreign[seg.source] = foreign.get(seg.source, 0) + hi - lo
    if host_bases:
        lo, hi, seg = max(host_bases, key=lambda t: t[1] - t[0])
        mapped_len = hi - lo
        if mapped_len >= cfg.min_mapped_fraction * length:
            host_start = seg.src_start + (lo - seg.start)
            left, right = lo - start, end - hi
            cigar = (f"{left}S" if left else "") + f"{mapped_len}M" + (f"{right}S" if right else "")
            return True, host_start, cigar, None
    taxid = max(foreign, key=foreign.get) if foreign else None
    if taxid and foreign[taxid] < cfg.min_mapped_fraction * length:
        taxid = None
    return False, -1, "*", taxid


class _Writer:
    def __init__(self, rng, cfg):
        self.rng, self.cfg = rng, cfg
        self.sam, self.kraken = [], []
        self.n = 0

    def _classify(self, name, taxid):
        if taxid is not None and self.rng.random() < self.cfg.kraken_sensitivity:
            self.kraken.append(f"C\t{name}\t{taxid}\t{self.cfg.read_length}\t{taxid}:1")
        else:
            self.kraken.append(f"U\t{name}\t0\t{self.cfg.read_length}\t0:1")

    def pair(self, chrom, r1, r2):
        """r1/r2: (mapped, host_start, cigar, taxid, seq, reverse)."""
        self.n += 1
        name = f"pair{self.n}"
        m1, m2 = r1[0], r2[0]
        mate_taxid = None
        for this, other, first in ((r1, r2, True), (r2, r1, False)):
            mapped, hstart, cigar, taxid, seq, rev = this
            flag = 0x1 | (0x40 if first else 0x80)
            if not mapped:
                flag |= 0x4
            if not other[0]:
                flag |= 0x8
            if rev:
                flag |= 0x10
            if other[5]:
                flag |= 0x20
            if mapped:
                rname, pos, mapq = chrom, hstart + 1, 60
            elif other[0]:
                rname, pos, mapq = chrom, other[1] + 1, 0
            else:
                rname, pos, mapq = "*", 0, 0
            qual = "I" * len(seq)
            self.sam.append(
                f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t=\t0\t0\t{seq}\t{qual}"
            )
            if not mapped:
                mate_taxid = taxid
        if not (m1 and m2):
            self._classify(name, mate_taxid)


def simulate(out_prefix, cfg=None):
    cfg = cfg or SimConfig()
    rng = random.Random(cfg.seed)
    donor = random_seq(rng, 200_000, gc=0.62)
    contaminant = random_seq(rng, 200_000, gc=0.33)
    writer = _Writer(rng, cfg)
    header = ["@HD\tVN:1.6\tSO:unsorted"]
    truth = []
    host_chroms = {}

    per_chrom = [cfg.n_insertions // cfg.n_chroms] * cfg.n_chroms
    for i in range(cfg.n_insertions % cfg.n_chroms):
        per_chrom[i] += 1

    L = cfg.read_length
    for c in range(cfg.n_chroms):
        chrom = f"chr{c + 1}"
        host = random_seq(rng, cfg.chrom_length, gc=0.41)
        host_chroms[chrom] = host
        header.append(f"@SQ\tSN:{chrom}\tLN:{len(host)}")
        sample, segments, chrom_truth = _build_sample_chrom(rng, host, donor, per_chrom[c], cfg, chrom)
        truth.extend(chrom_truth)
        n_pairs = int(cfg.coverage * len(sample) / (2 * L))
        for _ in range(n_pairs):
            frag = max(2 * L, int(rng.gauss(cfg.fragment_mean, cfg.fragment_sd)))
            start = rng.randint(0, len(sample) - frag)
            s1, s2 = start, start + frag - L
            p1 = _place(segments, s1, L, cfg)
            p2 = _place(segments, s2, L, cfg)
            seq1 = sample[s1 : s1 + L]
            seq2 = revcomp(sample[s2 : s2 + L])
            writer.pair(chrom, (*p1, seq1, False), (*p2, seq2, True))

    for _ in range(cfg.contaminant_pairs):
        start = rng.randint(0, len(contaminant) - 2 * L - 200)
        seq1 = contaminant[start : start + L]
        seq2 = revcomp(contaminant[start + 200 : start + 200 + L])
        cont = (False, -1, "*", str(CONTAMINANT_TAXID))
        writer.pair("*", (*cont, seq1, False), (*cont, seq2, True))

    chroms = list(host_chroms)
    for _ in range(cfg.chimera_pairs):
        chrom = rng.choice(chroms)
        host = host_chroms[chrom]
        hstart = rng.randint(0, len(host) - L)
        cstart = rng.randint(0, len(contaminant) - L)
        reverse = rng.random() < 0.5
        host_seq = host[hstart : hstart + L]
        r_host = (True, hstart, f"{L}M", None, revcomp(host_seq) if reverse else host_seq, reverse)
        r_cont = (False, -1, "*", str(CONTAMINANT_TAXID), contaminant[cstart : cstart + L], not reverse)
        writer.pair(chrom, r_host, r_cont)

    with open(f"{out_prefix}.sam", "w") as fh:
        fh.write("\n".join(header + writer.sam) + "\n")
    with open(f"{out_prefix}.kraken2.txt", "w") as fh:
        fh.write("\n".join(writer.kraken) + "\n")
    with open(f"{out_prefix}.truth.bed", "w") as fh:
        for chrom, pos, length in sorted(truth):
            fh.write(f"{chrom}\t{pos}\t{pos + 1}\t{DONOR_TAXID}\t{length}\n")
    return truth
