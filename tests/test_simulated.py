import csv

import pytest

from hgt_site.call import call_sites, write_table
from hgt_site.cli import main
from hgt_site.evaluate import evaluate, read_calls, read_truth
from hgt_site.simulate import CONTAMINANT_TAXID, HOST_TAXID, SimConfig, simulate


@pytest.fixture(scope="module")
def sim(tmp_path_factory):
    prefix = tmp_path_factory.mktemp("sim") / "sample"
    simulate(str(prefix), SimConfig(seed=7))
    return prefix


def test_recovers_true_sites(sim, tmp_path):
    out = tmp_path / "sites.tsv"
    sites = call_sites(f"{sim}.sam", f"{sim}.kraken2.txt", HOST_TAXID)
    write_table(sites, out)
    result = evaluate(read_calls(out), read_truth(f"{sim}.truth.bed"))
    assert result["recall"] == 1.0
    assert result["precision"] == 1.0


def test_chimeras_do_not_pass(sim, tmp_path):
    out = tmp_path / "sites.tsv"
    write_table(call_sites(f"{sim}.sam", f"{sim}.kraken2.txt", HOST_TAXID), out)
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    contaminant_rows = [r for r in rows if r["top_taxid"] == str(CONTAMINANT_TAXID)]
    assert contaminant_rows, "the simulation should produce some chimera clusters"
    assert all(r["filter"] != "PASS" for r in contaminant_rows)


def test_excluding_donor_removes_calls(sim):
    sites = call_sites(f"{sim}.sam", f"{sim}.kraken2.txt", HOST_TAXID, exclude=[562])
    assert all(s.stats["top_taxid"] != 562 for s in sites)


def test_cli_runs(sim, tmp_path):
    out = tmp_path / "cli.tsv"
    code = main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out),
    ])
    assert code == 0
    assert out.read_text().startswith("chrom\tstart")


def test_simulated_sam_has_mate_positions(sim):
    with open(f"{sim}.sam") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.split("\t")
            if not int(cols[1]) & 0x8:
                assert int(cols[7]) > 0, "a mapped mate must have a position"
