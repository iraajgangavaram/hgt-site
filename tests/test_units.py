import math

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
