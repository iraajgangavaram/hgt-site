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
