"""Classify unmapped reads against a donor sequence instead of a classifier.

Someone who already knows what the foreign DNA is, because they put it there,
should not have to build a Kraken2 database to find out where it went. Given a
FASTA of the donor, a vector, a transgene or a viral genome, every unmapped
read is checked against it directly and the ones that match are treated as
calls for the donor taxon.

Matching is exact k-mer matching in both orientations. For a donor of a few
kilobases to a few hundred kilobases the k-mer set is small enough to hold in
memory and fast enough to be beside the point. A 150 bp read holds 120 31-mers
and a sequencing error spoils at most 31 of them, so a handful of exact
matches is a low bar for a true donor read and a high one for anything else.

This is not an aligner and does not try to be. It answers "does this read come
from the donor", which is the only question the caller needs, and it will not
find a read that is diverged from the donor sequence throughout. For a donor
that differs from the reference by more than a few per cent, align the reads
and use a classifier.
"""

import gzip
import logging

# NCBI's taxon for a sequence that was built rather than sampled. It is the
# honest default for a vector or a construct, and --donor-taxid overrides it
# for anyone who knows the real taxon of their donor.
SYNTHETIC_CONSTRUCT = 32630

KMER = 31
MIN_HITS = 5

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

log = logging.getLogger(__name__)


def read_fasta(path):
    """Read a FASTA file into {name: sequence}, uppercased."""
    opener = gzip.open if str(path).endswith(".gz") else open
    sequences, name, chunks = {}, None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0] if len(line) > 1 else "", []
            else:
                chunks.append(line.upper())
    if name is not None:
        sequences[name] = "".join(chunks)
    if not sequences:
        raise ValueError(f"{path} holds no FASTA records.")
    return sequences


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def kmers(seq, k=KMER):
    """The k-mers of one sequence, skipping any that contain an N."""
    return {seq[i : i + k] for i in range(len(seq) - k + 1) if "N" not in seq[i : i + k]}


def build_index(sequences, k=KMER):
    """The k-mers of every donor record, on both strands.

    Both strands are stored rather than one canonical form of each k-mer. It
    doubles the size of the set, which is nothing for a donor of this size,
    and it means a read can be looked up as it is. Taking the reverse
    complement of every k-mer of every read instead was about three times
    slower, and reads are the side there are millions of.
    """
    index = set()
    for seq in sequences.values():
        index |= kmers(seq, k)
        index |= kmers(revcomp(seq), k)
    return index


def count_hits(seq, index, k=KMER, stop_at=None):
    """How many of a read's k-mers are in the donor.

    Counting stops once stop_at is reached, since the caller only needs to
    know whether the read clears the threshold.
    """
    hits = 0
    for i in range(len(seq) - k + 1):
        if seq[i : i + k] in index:
            hits += 1
            if stop_at is not None and hits >= stop_at:
                return hits
    return hits


def classify(read_seqs, index, taxid, k=KMER, min_hits=MIN_HITS):
    """Call the reads that match the donor. Returns {read name: taxid}."""
    calls = {}
    for name, seq in read_seqs.items():
        if not seq or len(seq) < k:
            continue
        if count_hits(seq.upper(), index, k, stop_at=min_hits) >= min_hits:
            calls[name] = taxid
    return calls


def classify_against_donor(path, read_seq_groups, taxid=SYNTHETIC_CONSTRUCT, k=KMER, min_hits=MIN_HITS):
    """Classify several groups of reads against a donor FASTA.

    read_seq_groups is a list of {read name: sequence} mappings, which are
    classified together so the donor is read and indexed once.
    """
    sequences = read_fasta(path)
    total = sum(len(s) for s in sequences.values())
    if total < k:
        raise ValueError(f"{path} is shorter than the {k} bp k-mer length.")
    index = build_index(sequences, k)
    log.info(
        "Donor %s: %d records, %d bp, %d distinct %d-mers.",
        path, len(sequences), total, len(index), k,
    )
    calls = {}
    for group in read_seq_groups:
        calls.update(classify(group, index, taxid, k, min_hits))
    return calls
