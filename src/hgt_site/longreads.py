"""Call sites from long reads, where the junction sits inside a single read.

Short reads come in pairs and hgt-site works from a pair whose mates land in
different places: one on the host, one foreign. A long read needs no partner.
A read that crosses an integration carries the host sequence and the inserted
sequence in one molecule, and an aligner places the host part and clips the
rest off, so the junction is the clip boundary and the clipped tail is the
foreign sequence itself.

That makes the long-read case the easier one in principle. There is no mate to
guess at, the junction is exact rather than an interval, and a chance chimera
has to join two molecules that happen to be long rather than two short ones,
which is rarer.

The evidence is reshaped into the same Anchor objects the paired path builds,
so clustering, the enrichment test, the flags and every output format are the
code that was already there rather than a second implementation of it. A read
clipped at both ends, which happens when a short element sits inside one read,
contributes an anchor for each end.

Classification is by `--donor` only. A classifier works on whole reads, and
here the foreign part is a tail of a read that aligned to the host, so a
per-read call would describe the host rather than the insertion. Matching the
clipped tail against the sequence you are looking for is the thing that
actually answers the question.
"""

from .alignments import Anchor
from .donor import KMER, MIN_HITS, build_index, count_hits, read_fasta

# A clip has to be long enough to carry real sequence rather than an adapter
# or a few mismatched bases at the end of an alignment.
MIN_CLIP = 200
# And the host side has to be substantial, or the read is not anchored
# anywhere in particular.
MIN_MAPPED = 500

SOFT_CLIP, HARD_CLIP = 4, 5


def clip_lengths(cigar):
    """Soft clip at the start and end of an alignment, in bases.

    Hard clips are reported as zero on purpose. The sequence of a hard clip is
    not in the record, so there is nothing to match against the donor, and a
    supplementary alignment carrying it is skipped anyway.
    """
    if not cigar:
        return 0, 0
    first, last = cigar[0], cigar[-1]
    start = first[1] if first[0] == SOFT_CLIP else 0
    end = last[1] if last[0] == SOFT_CLIP and len(cigar) > 1 else 0
    return start, end


def _clips_of(read, min_clip, min_mapped):
    """The clipped tails of one read that are long enough to be worth checking.

    Yields (side, sequence). `side` is "left" when the clip is before the
    aligned part, which puts the junction at the alignment's start.
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    if not cigar or not seq:
        return
    if read.reference_length is None or read.reference_length < min_mapped:
        return
    start_clip, end_clip = clip_lengths(cigar)
    if start_clip >= min_clip:
        yield "left", seq[:start_clip]
    if end_clip >= min_clip:
        yield "right", seq[len(seq) - end_clip :]


def collect_clipped(
    path, donor, taxid, min_mapq=20, min_clip=MIN_CLIP, min_mapped=MIN_MAPPED,
    reference=None, threads=1, index=None, k=KMER, min_hits=MIN_HITS,
):
    """Anchors from long reads clipped onto sequence that matches the donor.

    Returns the anchors, the taxid call for each one, the clipped sequence for
    each one, and the contig lengths, which is what the paired path hands on
    to clustering and scoring.
    """
    import pysam

    donor_index = build_index(read_fasta(donor), k)

    mode = "rc" if str(path).endswith(".cram") else "rb"
    anchors, calls, seqs = [], {}, {}
    reads = clipped = foreign = 0
    with pysam.AlignmentFile(
        str(path), mode, reference_filename=reference, threads=threads,
        index_filename=index,
    ) as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            reads += 1
            if read.mapping_quality < min_mapq:
                continue
            for side, tail in _clips_of(read, min_clip, min_mapped):
                clipped += 1
                if count_hits(tail, donor_index, k, stop_at=min_hits) < min_hits:
                    continue
                foreign += 1
                # One read can cross two junctions, so each end gets its own
                # name. Everything downstream keys on the anchor name.
                name = f"{read.query_name}:{side[0]}"
                # The anchor is the junction, not the read. A paired anchor
                # spans a read length and the junction is somewhere past its
                # end, so the site inherits that uncertainty. Here the clip
                # boundary is the junction, to the base. Carrying the whole
                # alignment instead would make a site as wide as the reads
                # that built it, tens of kilobases, and the enrichment test
                # divides by the genome in proportion to that width, so wide
                # sites are tested against a far larger expected count and
                # real ones fail.
                junction = (
                    read.reference_start if side == "left" else read.reference_end - 1
                )
                anchors.append(
                    Anchor(
                        name=name,
                        chrom=read.reference_name,
                        start=junction,
                        end=junction + 1,
                        # A left clip puts the junction at the start of the
                        # alignment, which is what a reverse anchor means to
                        # the clustering code.
                        reverse=(side == "left"),
                        mapq=read.mapping_quality,
                    )
                )
                calls[name] = taxid
                seqs[name] = tail

    stats = {
        "long_reads": reads,
        "clipped": clipped,
        "foreign_clips": foreign,
    }
    return anchors, calls, seqs, lengths, stats
