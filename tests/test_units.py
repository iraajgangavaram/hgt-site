import math

import pytest

from hgt_site.alignments import Anchor
from hgt_site.kraken import TaxonFilter, normalise_read_name, parse_taxid
from hgt_site.sites import benjamini_hochberg, cluster_anchors, poisson_sf, sequence_complexity


def test_read_name_suffix():
    assert normalise_read_name("r1/1") == "r1"
    assert normalise_read_name("r1/2") == "r1"
    assert normalise_read_name("r1") == "r1"


def test_parse_taxid_with_names():
    assert parse_taxid("562") == 562
    assert parse_taxid("Escherichia coli (taxid 562)") == 562


def test_taxon_filter_with_lineage(tmp_path):
    nodes = tmp_path / "nodes.dmp"
    # 1 root > 131567 > 2759 Eukaryota > 9606 human; 131567 > 2 Bacteria > 562
    rows = [(1, 1), (131567, 1), (2759, 131567), (9606, 2759), (2, 131567), (562, 2)]
    nodes.write_text("".join(f"{t}\t|\t{p}\t|\tno rank\t|\n" for t, p in rows))
    f = TaxonFilter(9606, nodes)
    assert f.is_foreign(562)
    assert f.is_foreign(2)
    assert not f.is_foreign(9606)
    assert not f.is_foreign(2759)  # ancestor of the host, uninformative
    assert not f.is_foreign(1)


def test_excluded_taxid_covers_its_clade(tmp_path):
    """Excluding a family excludes the species in it, once a taxonomy is given.

    Matching only the exact taxid made the option useless for the case it is
    needed in: a bacterial host with relatives in the sample, where the reads
    are called at species level and the thing you want to exclude is the
    family.
    """
    nodes = tmp_path / "nodes.dmp"
    # 2 Bacteria > 1224 Proteobacteria > 91347 Enterobacterales >
    #   543 Enterobacteriaceae > {562 E. coli, 28901 Salmonella, 573 Klebsiella}
    # and 1224 > 72274 Pseudomonadales > 286 Pseudomonas
    rows = [
        (1, 1), (2, 1), (1224, 2), (91347, 1224), (543, 91347),
        (562, 543), (28901, 543), (573, 543),
        (72274, 1224), (286, 72274),
    ]
    nodes.write_text("".join(f"{a}\t|\t{b}\t|\tno rank\t|\n" for a, b in rows))

    plain = TaxonFilter(562, nodes)
    assert plain.is_foreign(28901) and plain.is_foreign(573)

    excluded = TaxonFilter(562, nodes, exclude=[543])
    assert not excluded.is_foreign(28901), "Salmonella is inside the excluded family"
    assert not excluded.is_foreign(573), "Klebsiella is inside the excluded family"
    assert not excluded.is_foreign(543)
    assert excluded.is_foreign(286), "Pseudomonas is outside it and must survive"


def test_excluded_taxid_is_exact_without_a_taxonomy():
    """With no nodes.dmp there is no clade to expand, so it matches exactly."""
    f = TaxonFilter(562, exclude=[543])
    assert not f.is_foreign(543)
    assert f.is_foreign(28901)


def test_clip_lengths_reads_both_ends_and_ignores_hard_clips():
    """Soft clips carry sequence to match against the donor; hard clips do not."""
    from hgt_site.longreads import clip_lengths

    MATCH, SOFT, HARD = 0, 4, 5
    assert clip_lengths([(SOFT, 300), (MATCH, 5000)]) == (300, 0)
    assert clip_lengths([(MATCH, 5000), (SOFT, 420)]) == (0, 420)
    assert clip_lengths([(SOFT, 300), (MATCH, 5000), (SOFT, 420)]) == (300, 420)
    assert clip_lengths([(MATCH, 5000)]) == (0, 0)
    assert clip_lengths([]) == (0, 0)
    # A hard clip has no sequence in the record, so it is worth nothing here.
    assert clip_lengths([(HARD, 900), (MATCH, 5000), (HARD, 900)]) == (0, 0)
    # A single operation must not be counted at both ends.
    assert clip_lengths([(SOFT, 300)]) == (300, 0)


def test_long_read_clips_are_taken_from_the_right_end():
    """The left clip is the head of the read and the right clip is the tail."""
    from hgt_site.longreads import _clips_of

    class FakeRead:
        cigartuples = [(4, 10), (0, 1000), (4, 20)]
        query_sequence = "L" * 10 + "M" * 1000 + "R" * 20
        reference_length = 1000

    clips = dict(_clips_of(FakeRead(), min_clip=5, min_mapped=100))
    assert clips["left"] == "L" * 10
    assert clips["right"] == "R" * 20
    # Too short a host side means the read is not anchored anywhere.
    assert dict(_clips_of(FakeRead(), min_clip=5, min_mapped=5000)) == {}
    # Too short a clip is not worth matching against the donor.
    assert dict(_clips_of(FakeRead(), min_clip=50, min_mapped=100)) == {}


def test_control_marks_shared_positions_only(tmp_path):
    """A site the control also has is flagged; one it does not have is left alone."""
    from hgt_site.control import flag_in_control, in_control, load_control
    from hgt_site.sites import Site

    table = tmp_path / "normal.tsv"
    table.write_text(
        "chrom\tstart\tend\tn_pairs\tfilter\n"
        # Enough reads to be reproducible, and it failed its own filters,
        # which is still evidence about the tumour.
        "chr1\t1000\t1200\t7\tLOW_SUPPORT\n"
        "chr1\t5000\t5200\t9\tPASS\n"
        # A couple of reads is what chance chimeras give everywhere, so these
        # must not count or the control rejects the whole genome.
        "chr1\t8000\t8200\t2\tLOW_SUPPORT\n"
        "chr2\t9000\t9200\t0\tLOW_SUPPORT\n"
    )
    positions = load_control(table)
    assert positions == {"chr1": [1100, 5100]}, "weak candidates are noise, not evidence"

    assert in_control(positions, "chr1", 1050, 1150)
    assert in_control(positions, "chr1", 5100, 5300)
    # Just inside and just outside the 300 bp window.
    assert in_control(positions, "chr1", 1350, 1380)
    assert not in_control(positions, "chr1", 1450, 1500)
    assert not in_control(positions, "chr1", 8000, 8200), "two reads must not count"
    assert not in_control(positions, "chr3", 1100, 1200)

    sites = [Site("chr1", 1050, 1150), Site("chr1", 40000, 40100)]
    assert flag_in_control(sites, positions) == 1
    assert sites[0].flags == ["IN_CONTROL"]
    assert sites[1].flags == []


def test_control_rejects_a_file_that_is_not_a_table(tmp_path):
    from hgt_site.control import load_control

    bad = tmp_path / "notes.txt"
    bad.write_text("this is not a sites table\n")
    with pytest.raises(ValueError, match="hgt-site table"):
        load_control(bad)


def test_complexity():
    assert sequence_complexity("A" * 100) < 0.05
    assert sequence_complexity("ACGTTGCAAGCTTCGAGATCCGATAGCTAGCATGCATCGA") > 0.8


def test_poisson_sf():
    assert poisson_sf(0, 3.0) == 1.0
    assert math.isclose(poisson_sf(1, 2.0), 1 - math.exp(-2.0))


def test_bh_monotone():
    q = benjamini_hochberg([0.01, 0.04, 0.03, 0.5])
    assert q[0] <= q[2] <= q[1] <= q[3]


def test_cluster_uses_junction_side():
    anchors = [
        Anchor("a", "chr1", 1000, 1100, False, 60),  # junction side 1100
        Anchor("b", "chr1", 1400, 1500, True, 60),  # junction side 1400
        Anchor("c", "chr1", 9000, 9100, False, 60),
    ]
    sites = cluster_anchors(anchors, max_gap=600, max_span=1000)
    assert [len(s.anchors) for s in sites] == [2, 1]


def test_cluster_span_limit_stops_chaining():
    anchors = [Anchor(str(i), "chr1", p, p + 100, False, 60) for i, p in enumerate(range(0, 5000, 500))]
    sites = cluster_anchors(anchors, max_gap=600, max_span=1000)
    assert max(len(s.anchors) for s in sites) <= 3


def test_mixed_taxa_allows_background_chimeras():
    from collections import Counter

    from hgt_site.call import DEFAULTS
    from hgt_site.sites import score_sites

    # 14 donor pairs plus 4 random chimeras of another taxon in one window.
    anchors = [Anchor(f"d{i}", "chr1", 5000 + i * 10, 5100 + i * 10, i % 2 == 0, 60) for i in range(14)]
    anchors += [Anchor(f"c{i}", "chr1", 5050 + i * 20, 5150 + i * 20, False, 60) for i in range(4)]
    calls = {f"d{i}": 562 for i in range(14)}
    calls.update({f"c{i}": 573 for i in range(4)})
    # Genome-wide, the other taxon has many scattered chimeric anchors.
    scattered = [Anchor(f"x{i}", "chr1", 20000 + i * 300, 20100 + i * 300, False, 60) for i in range(30000)]
    calls.update({f"x{i}": 573 for i in range(30000)})
    sites = cluster_anchors(anchors + scattered, max_gap=300, max_span=800)
    sites = score_sites(sites, calls, {}, Counter(), 10_000_000, 600, DEFAULTS)
    real = max(sites, key=lambda s: len(s.anchors))
    assert real.stats["top_taxid_fraction"] < 0.8
    assert "MIXED_TAXA" not in real.flags


def test_mixed_taxa_flags_excess_other_taxa():
    from collections import Counter

    from hgt_site.call import DEFAULTS
    from hgt_site.sites import score_sites

    # Same site, but the other taxon is rare genome-wide, so 4 reads here is an excess.
    anchors = [Anchor(f"d{i}", "chr1", 5000 + i * 10, 5100 + i * 10, i % 2 == 0, 60) for i in range(14)]
    anchors += [Anchor(f"c{i}", "chr1", 5050 + i * 20, 5150 + i * 20, False, 60) for i in range(4)]
    calls = {f"d{i}": 562 for i in range(14)}
    calls.update({f"c{i}": 573 for i in range(4)})
    scattered = [Anchor(f"x{i}", "chr1", 20000 + i * 3000, 20100 + i * 3000, False, 60) for i in range(3000)]
    calls.update({f"x{i}": 573 for i in range(3000)})
    sites = cluster_anchors(anchors + scattered, max_gap=300, max_span=800)
    sites = score_sites(sites, calls, {}, Counter(), 10_000_000, 600, DEFAULTS)
    real = max(sites, key=lambda s: len(s.anchors))
    assert "MIXED_TAXA" in real.flags


def test_only_needed_records_are_kept():
    from hgt_site.alignments import _wanted

    assert _wanted(0x1 | 0x8)  # mapped read, mate unmapped
    assert _wanted(0x1 | 0x4)  # unmapped read, mate mapped
    assert _wanted(0x1 | 0x4 | 0x8)  # both unmapped
    assert not _wanted(0x1)  # both mates mapped
    assert not _wanted(0x8)  # not paired
    assert not _wanted(0x1 | 0x8 | 0x100)  # secondary
    assert not _wanted(0x1 | 0x8 | 0x800)  # supplementary
    assert not _wanted(0x1 | 0x8 | 0x400)  # duplicate


def test_limits_widen_only_for_long_fragments():
    from hgt_site.alignments import ReadPairs
    from hgt_site.call import DEFAULTS, fit_limits

    def pairs(fragment):
        return ReadPairs([], {}, set(), {}, [fragment] * 2000, [100] * 2000)

    assert fit_limits(pairs(350), dict(DEFAULTS), set()) == DEFAULTS
    wide = fit_limits(pairs(900), dict(DEFAULTS), set())
    assert wide["max_span"] == 2 * (900 - 100)
    assert wide["max_gap"] > DEFAULTS["max_gap"] and wide["window"] > DEFAULTS["window"]
    # Limits the user set are left alone.
    assert fit_limits(pairs(900), dict(DEFAULTS), {"max_span"}) == DEFAULTS
    # Too few pairs to estimate: defaults are kept.
    assert fit_limits(ReadPairs([], {}, set(), {}, [900] * 10, [100] * 10), dict(DEFAULTS), set()) == DEFAULTS
