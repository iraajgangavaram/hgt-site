import csv
import random

import pysam

from hgt_site.longreads import collect_clipped
from hgt_site.call import call_sites, write_table


def _read(name, sequence, start, cigar, mapq=60):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = 0
    read.reference_id = 0
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigartuples = cigar
    read.query_sequence = sequence
    read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    return read


def test_collect_clipped_matches_donor_tails_at_both_junction_sides(tmp_path):
    rng = random.Random(41)
    donor_sequence = "".join(rng.choices("ACGT", k=1_000))
    donor = tmp_path / "donor.fa"
    donor.write_text(f">donor\n{donor_sequence}\n")
    host_sequence = "".join(rng.choices("ACGT", k=1_200))
    alignment = tmp_path / "reads.bam"

    with pysam.AlignmentFile(
        alignment,
        "wb",
        header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 5_000}]},
    ) as bam:
        bam.write(_read("right", host_sequence[:600] + donor_sequence[:250], 1_000, [(0, 600), (4, 250)]))
        bam.write(_read("left", donor_sequence[300:550] + host_sequence[600:], 2_000, [(4, 250), (0, 600)]))
        bam.write(_read("noise", host_sequence[:600] + "T" * 250, 3_000, [(0, 600), (4, 250)]))

    anchors, calls, sequences, lengths, stats = collect_clipped(alignment, donor, 562)

    assert [(a.name, a.start, a.end, a.reverse) for a in anchors] == [
        ("right:r", 1_599, 1_600, False),
        ("left:l", 2_000, 2_001, True),
    ]
    assert calls == {"right:r": 562, "left:l": 562}
    assert sequences == {"right:r": donor_sequence[:250], "left:l": donor_sequence[300:550]}
    assert lengths == {"chr1": 5_000}
    assert stats == {"long_reads": 3, "clipped": 3, "foreign_clips": 2}


def test_long_read_call_reports_supported_donor_junction(tmp_path):
    rng = random.Random(93)
    donor_sequence = "".join(rng.choices("ACGT", k=1_000))
    donor = tmp_path / "donor.fa"
    donor.write_text(f">donor\n{donor_sequence}\n")
    host_sequence = "".join(rng.choices("ACGT", k=600))
    alignment = tmp_path / "reads.bam"

    with pysam.AlignmentFile(
        alignment,
        "wb",
        header={"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 5_000}]},
    ) as bam:
        for i in range(8):
            bam.write(
                _read(
                    f"junction{i}",
                    host_sequence + donor_sequence[:250],
                    1_000,
                    [(0, 600), (4, 250)],
                )
            )

    sites = call_sites(
        alignment,
        None,
        host_taxid=9606,
        donor=donor,
        long_reads=True,
        split_reads=False,
    )

    assert len(sites) == 1
    assert sites[0].start == 1_599
    assert sites[0].end == 1_600
    output = tmp_path / "sites.tsv"
    write_table(sites, output)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["filter"] == "PASS"
    assert rows[0]["flags"] == "ONE_SIDED"
