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


@pytest.fixture(scope="module")
def sim_split(tmp_path_factory):
    # Reads may align with as little as 30% of their length, so many reads
    # crossing a junction are kept as soft-clipped split reads.
    prefix = tmp_path_factory.mktemp("split") / "sample"
    simulate(str(prefix), SimConfig(seed=7, coverage=20, min_mapped_fraction=0.3))
    return prefix


def _breakpoints_match_truth(sites, truth_bed):
    truth = {}
    for line in open(truth_bed):
        chrom, start, *_ = line.split("\t")
        truth.setdefault(chrom, []).append(int(start))
    passed = [s for s in sites if not [f for f in s.flags if f != "ONE_SIDED"]]
    assert len(passed) == 5
    for site in passed:
        st = site.stats
        assert st["left_split_reads"] >= 2 or st["right_split_reads"] >= 2
        for key in ("left_breakpoint", "right_breakpoint"):
            if st[key] != "":
                assert st[key] in truth[site.chrom], (key, st[key], truth[site.chrom])


def test_split_reads_give_exact_junctions_sam(sim_split):
    sites = call_sites(f"{sim_split}.sam", f"{sim_split}.kraken2.txt", HOST_TAXID)
    _breakpoints_match_truth(sites, f"{sim_split}.truth.bed")


def _rows(sites):
    """Sites as comparable rows, so two formats can be compared whatever the order."""
    return sorted(
        (s.chrom, s.start, s.end, tuple(sorted((k, str(v)) for k, v in s.stats.items())))
        for s in sites
    )


def _same_sites(a, b, label_a, label_b):
    rows_a, rows_b = _rows(a), _rows(b)
    only_a = [r for r in rows_a if r not in rows_b]
    only_b = [r for r in rows_b if r not in rows_a]
    assert not only_a and not only_b, (
        f"{len(only_a)} sites only from {label_a}, {len(only_b)} only from {label_b}; "
        f"first from {label_a}: {only_a[:1]}; first from {label_b}: {only_b[:1]}"
    )


def test_split_reads_give_exact_junctions_bam(sim_split, tmp_path):
    pysam = pytest.importorskip("pysam")
    bam = tmp_path / "sample.bam"
    pysam.sort("-o", str(bam), f"{sim_split}.sam")
    pysam.index(str(bam))
    sites = call_sites(str(bam), f"{sim_split}.kraken2.txt", HOST_TAXID)
    _breakpoints_match_truth(sites, f"{sim_split}.truth.bed")
    sam_sites = call_sites(f"{sim_split}.sam", f"{sim_split}.kraken2.txt", HOST_TAXID)
    _same_sites(sites, sam_sites, "the BAM file", "the SAM file")


def test_unindexed_bam_falls_back_to_pairs(sim_split, tmp_path, caplog):
    pysam = pytest.importorskip("pysam")
    bam = tmp_path / "unsorted.bam"
    pysam.view("-b", "-o", str(bam), f"{sim_split}.sam", catch_stdout=False)
    with caplog.at_level("INFO", logger="hgt_site"):
        sites = call_sites(str(bam), f"{sim_split}.kraken2.txt", HOST_TAXID)
    assert "no index" in caplog.text
    assert all(s.stats["left_split_reads"] == 0 and s.stats["right_split_reads"] == 0 for s in sites)
    assert sum(1 for s in sites if not [f for f in s.flags if f != "ONE_SIDED"]) == 5


def test_split_reads_can_be_turned_off(sim_split):
    sites = call_sites(f"{sim_split}.sam", f"{sim_split}.kraken2.txt", HOST_TAXID, split_reads=False)
    assert all(s.stats["left_breakpoint"] == "" for s in sites)


def test_vcf_uses_exact_junction(sim_split, tmp_path):
    vcf = tmp_path / "s.vcf"
    assert main([
        "call", "--bam", f"{sim_split}.sam", "--kraken", f"{sim_split}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(tmp_path / "s.tsv"), "--vcf", str(vcf),
    ]) == 0
    truth = {int(line.split("\t")[1]) for line in open(f"{sim_split}.truth.bed")}
    passed = [line.split("\t") for line in vcf.read_text().splitlines() if not line.startswith("#")]
    passed = [r for r in passed if r[6] == "PASS"]
    assert len(passed) == 5
    for r in passed:
        assert "IMPRECISE" not in r[7]
        assert int(r[1]) in truth  # POS is the last host base before the insertion


def test_multiqc_table(sim, tmp_path):
    import json

    mqc = tmp_path / "sample_mqc.json"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(tmp_path / "s.tsv"), "--multiqc", str(mqc),
    ]) == 0
    content = json.loads(mqc.read_text())
    assert content["plot_type"] == "table"
    row = content["data"]["sample"]
    assert row["passed_sites"] == 5
    assert row["candidate_sites"] >= row["passed_sites"]
    assert row["foreign_anchored_pairs"] <= row["anchored_pairs_with_call"] <= row["anchored_pairs"]
    assert set(content["headers"]) <= set(row)


def test_cram_input(sim_split, tmp_path):
    pysam = pytest.importorskip("pysam")
    ref = f"{sim_split}.host.fa"
    bam, cram = tmp_path / "s.bam", tmp_path / "s.cram"
    pysam.sort("-o", str(bam), f"{sim_split}.sam")
    pysam.view("-C", "-T", ref, "-o", str(cram), str(bam), catch_stdout=False)
    pysam.index(str(cram))
    from_cram = call_sites(str(cram), f"{sim_split}.kraken2.txt", HOST_TAXID, reference=ref)
    from_sam = call_sites(f"{sim_split}.sam", f"{sim_split}.kraken2.txt", HOST_TAXID)
    _same_sites(from_cram, from_sam, "the CRAM file", "the SAM file")


def test_empty_kraken_file(sim, tmp_path, caplog):
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    with caplog.at_level("WARNING", logger="hgt_site"):
        assert call_sites(f"{sim}.sam", str(empty), HOST_TAXID) == []
    assert "no classified reads" in caplog.text


def test_outputs_are_valid_when_nothing_is_found(sim, tmp_path):
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    out, bed, vcf = tmp_path / "s.tsv", tmp_path / "s.bed", tmp_path / "s.vcf"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", str(empty), "--host-taxid", str(HOST_TAXID),
        "-o", str(out), "--bed", str(bed), "--vcf", str(vcf),
    ]) == 0
    assert out.read_text().count("\n") == 1  # header only
    assert bed.read_text() == ""
    assert vcf.read_text().splitlines()[-1].startswith("#CHROM")


def test_results_do_not_depend_on_read_order(sim_split, tmp_path):
    """The same reads in a different order must give the same sites."""
    import random

    lines = open(f"{sim_split}.sam").read().splitlines(keepends=True)
    header = [line for line in lines if line.startswith("@")]
    body = [line for line in lines if not line.startswith("@")]
    random.Random(3).shuffle(body)
    shuffled = tmp_path / "shuffled.sam"
    shuffled.write_text("".join(header + body))
    a = call_sites(f"{sim_split}.sam", f"{sim_split}.kraken2.txt", HOST_TAXID)
    b = call_sites(str(shuffled), f"{sim_split}.kraken2.txt", HOST_TAXID)
    assert [(s.chrom, s.start, s.end, s.stats) for s in a] == [(s.chrom, s.start, s.end, s.stats) for s in b]
    _same_sites(a, b, "the original order", "the shuffled order")
