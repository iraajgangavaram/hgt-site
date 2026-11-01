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


def test_warns_when_read_names_do_not_match(sim, tmp_path, caplog):
    renamed = tmp_path / "renamed.kraken2.txt"
    with open(f"{sim}.kraken2.txt") as src, open(renamed, "w") as dst:
        for line in src:
            cols = line.split("\t")
            cols[1] = "other_" + cols[1]
            dst.write("\t".join(cols))
    with caplog.at_level("WARNING", logger="hgt_site"):
        sites = call_sites(f"{sim}.sam", str(renamed), HOST_TAXID)
    assert sites == []
    assert "read names match" in caplog.text


def test_bed_has_one_line_per_passed_site(sim, tmp_path):
    out, bed = tmp_path / "s.tsv", tmp_path / "s.bed"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out), "--bed", str(bed),
    ]) == 0
    with open(out) as handle:
        passed = [r for r in csv.DictReader(handle, delimiter="\t") if r["filter"] == "PASS"]
    lines = bed.read_text().splitlines()
    assert len(lines) == len(passed) == 5
    for line in lines:
        chrom, start, end, name, score, strand = line.split("\t")
        assert int(start) < int(end) and name.startswith("taxid_") and strand == "."


def test_kraken_calls_can_be_split_across_files(sim, tmp_path):
    lines = open(f"{sim}.kraken2.txt").read().splitlines(keepends=True)
    half = len(lines) // 2
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    a.write_text("".join(lines[:half]))
    b.write_text("".join(lines[half:]))
    whole = call_sites(f"{sim}.sam", f"{sim}.kraken2.txt", HOST_TAXID)
    split = call_sites(f"{sim}.sam", [str(a), str(b)], HOST_TAXID)
    assert [s.stats for s in whole] == [s.stats for s in split]


def test_vcf_output(sim, tmp_path):
    out, vcf = tmp_path / "s.tsv", tmp_path / "s.vcf"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out), "--vcf", str(vcf),
    ]) == 0
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    lines = vcf.read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    records = [line.split("\t") for line in lines if not line.startswith("#")]
    assert len(records) == len(rows)
    assert sum(r[6] == "PASS" for r in records) == sum(r["filter"] == "PASS" for r in rows) == 5
    positions = [(r[0], int(r[1])) for r in records]
    assert positions == sorted(positions)
    for r in records:
        assert r[3] == "N" and r[4] == "<INS>" and "SVTYPE=INS" in r[7]
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(vcf)) as handle:
        parsed = list(handle)
    assert len(parsed) == len(rows)
