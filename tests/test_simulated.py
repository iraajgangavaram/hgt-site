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


def test_vcf_header_records_how_it_was_made(sim, tmp_path):
    vcf = tmp_path / "s.vcf"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(tmp_path / "s.tsv"), "--vcf", str(vcf),
        "--sample-name", "testsample", "--min-support", "4",
    ]) == 0
    text = vcf.read_text()
    # Contig lengths come from the alignment header, so other tools can index.
    assert "##contig=<ID=chr1,length=" in text
    assert "##hgt_siteCommand=hgt-site call " in text
    # The limits actually used, not the defaults.
    assert "##hgt_siteSetting=<ID=min_support,Value=4>" in text
    header = next(line for line in text.splitlines() if line.startswith("#CHROM"))
    assert header.split("\t")[8:] == ["FORMAT", "testsample"]


def test_vcf_quality_and_confidence_interval(sim, tmp_path):
    vcf = tmp_path / "s.vcf"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(tmp_path / "s.tsv"), "--vcf", str(vcf),
    ]) == 0
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(vcf)) as handle:
        records = list(handle)
    passed = [r for r in records if list(r.filter.keys()) == ["PASS"]]
    assert len(passed) == 5
    for record in passed:
        # A site that passed the enrichment test is well above the q-value
        # threshold of 0.05, which is a QUAL of about 13.
        assert record.qual > 13
        assert record.samples[0]["PR"] == record.info["NPAIRS"]
        assert record.samples[0]["GT"] == (None, None)
    failed = [r for r in records if "NOT_ENRICHED" in r.filter.keys()]
    assert failed and all(r.qual is not None for r in failed)
    # Almost every junction here is estimated from read pairs rather than from
    # split reads, and an estimated one has to carry an interval that contains
    # the position reported.
    imprecise = [r for r in records if "IMPRECISE" in r.info]
    assert len(imprecise) > 0.9 * len(records)
    for record in imprecise:
        low, high = record.info["CIPOS"]
        assert low <= 0 <= high
        assert record.info["CIPOS"] == record.info["CIEND"]
    # A junction found by split reads is exact, so it carries no interval.
    for record in records:
        if "IMPRECISE" not in record.info:
            assert "CIPOS" not in record.info
            assert record.samples[0]["SR"] > 0


def test_vcf_gz_is_compressed_and_indexed(sim, tmp_path):
    pysam = pytest.importorskip("pysam")
    vcf = tmp_path / "s.vcf.gz"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(tmp_path / "s.tsv"), "--vcf", str(vcf),
    ]) == 0
    assert vcf.exists() and (tmp_path / "s.vcf.gz.tbi").exists()
    assert not (tmp_path / "s.vcf").exists()
    assert vcf.read_bytes()[:2] == b"\x1f\x8b"
    with pysam.VariantFile(str(vcf)) as handle:
        # A region query only works on an indexed file.
        assert list(handle.fetch("chr1"))


def _truth_positions(prefix):
    rows = [line.split("\t") for line in open(f"{prefix}.truth.bed")]
    return [(r[0], int(r[1])) for r in rows]


def test_annotation_names_the_gene(sim, tmp_path):
    chrom, position = _truth_positions(sim)[0]
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "#!comment\n"
        # A gene containing the first real insertion.
        f'{chrom}\tsrc\tgene\t{position - 800}\t{position + 800}\t.\t+\t.\t'
        'gene_id "G1"; gene_name "INSIDE";\n'
        # A gene well downstream of it, and a transcript line that must be ignored.
        f'{chrom}\tsrc\tgene\t{position + 5000}\t{position + 6000}\t.\t-\t.\t'
        'gene_id "G2"; gene_name "NEARBY";\n'
        f'{chrom}\tsrc\ttranscript\t{position - 800}\t{position + 800}\t.\t+\t.\t'
        'gene_id "G1"; transcript_id "T1";\n'
    )
    out = tmp_path / "s.tsv"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out), "--annotation", str(gtf),
    ]) == 0
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    inside = [r for r in rows if r["gene"] == "INSIDE"]
    assert inside, "the site at a known insertion should be inside the gene covering it"
    assert all(r["gene_id"] == "G1" for r in inside)
    assert any(r["gene_distance"] == "0" for r in inside)
    # Every row has the columns, whether or not a gene was found.
    assert all("gene" in r and "gene_distance" in r for r in rows)
    # A site can sit inside the second gene too, so the only thing that has to
    # hold is that a distance is reported whenever a gene is.
    named = [r for r in rows if r["gene"]]
    assert named and all(r["gene_distance"] != "" for r in named)


def test_annotation_respects_the_distance_limit(sim, tmp_path):
    chrom, position = _truth_positions(sim)[0]
    gtf = tmp_path / "far.gtf"
    gtf.write_text(
        f'{chrom}\tsrc\tgene\t{position + 50000}\t{position + 51000}\t.\t+\t.\t'
        'gene_id "G9"; gene_name "FAR";\n'
    )
    out = tmp_path / "s.tsv"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out),
        "--annotation", str(gtf), "--max-gene-distance", "1000",
    ]) == 0
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    named = [r for r in rows if r["gene"]]
    assert all(int(r["gene_distance"]) <= 1000 for r in named)
    # The known insertion is 50 kb from the only gene, so it gets none.
    at_truth = [r for r in rows if int(r["start"]) <= position <= int(r["end"])]
    assert at_truth and all(r["gene"] == "" for r in at_truth)


def test_annotation_reads_gff3_and_rejects_a_file_with_no_genes(sim, tmp_path):
    from hgt_site.annotate import load_genes

    chrom, position = _truth_positions(sim)[0]
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "##gff-version 3\n"
        f"{chrom}\tEnsembl\tgene\t{position - 100}\t{position + 100}\t.\t+\t.\t"
        "ID=gene:ENSG00000000001;Name=TP53LIKE\n"
    )
    genes = load_genes(str(gff))
    starts, records = genes[chrom]
    assert records[0][2] == "TP53LIKE"
    assert records[0][3] == "ENSG00000000001"  # the "gene:" prefix is stripped

    empty = tmp_path / "none.gtf"
    empty.write_text("#!only a comment\n")
    with pytest.raises(ValueError, match="No gene records"):
        load_genes(str(empty))


def _as_format(kraken_path, out_path, fmt):
    """Rewrite Kraken2 per-read output in another classifier's layout."""
    rows = []
    for line in open(kraken_path):
        cols = line.rstrip("\n").split("\t")
        if cols and cols[0] == "C":
            rows.append((cols[1], cols[2]))
    with open(out_path, "w") as handle:
        if fmt == "centrifuge":
            handle.write("readID\tseqID\ttaxID\tscore\t2ndBestScore\thitLength\tqueryLength\tnumMatches\n")
        for name, taxid in rows:
            if fmt == "kaiju":
                handle.write(f"C\t{name}\t{taxid}\n")
            elif fmt == "centrifuge":
                handle.write(f"{name}\tseq1\t{taxid}\t400\t0\t150\t150\t1\n")
            elif fmt == "diamond":
                handle.write(f"{name}\t{taxid}\t1e-40\n")
            elif fmt == "blast":
                handle.write(f"{name}\tref1\t98.5\t150\t2\t0\t1\t150\t1\t150\t1e-60\t250\t{taxid}\n")
            else:
                handle.write(f"{name}\t{taxid}\n")
    return out_path


@pytest.mark.parametrize("fmt", ["kaiju", "centrifuge", "diamond", "blast", "tsv"])
def test_other_classifiers_give_the_same_sites(sim, tmp_path, fmt):
    from hgt_site.classifiers import detect_format

    converted = _as_format(f"{sim}.kraken2.txt", tmp_path / f"calls.{fmt}", fmt)
    expected = "kraken2" if fmt == "kaiju" else fmt
    assert detect_format(converted) == expected

    baseline, other = tmp_path / "kraken.tsv", tmp_path / "other.tsv"
    for output, calls in ((baseline, f"{sim}.kraken2.txt"), (other, str(converted))):
        assert main([
            "call", "--bam", f"{sim}.sam", "--classifications", str(calls),
            "--host-taxid", str(HOST_TAXID), "-o", str(output),
        ]) == 0
    assert other.read_text() == baseline.read_text()


def test_classifier_can_be_given_explicitly_and_reads_gzip(sim, tmp_path):
    import gzip

    from hgt_site.classifiers import load_classifications

    plain = _as_format(f"{sim}.kraken2.txt", tmp_path / "calls.tsv", "tsv")
    packed = tmp_path / "calls.tsv.gz"
    packed.write_bytes(gzip.compress(plain.read_bytes()))
    assert load_classifications([str(packed)], "tsv") == load_classifications([str(plain)])


def test_unreadable_classifier_output_says_so(tmp_path):
    from hgt_site.classifiers import detect_format

    odd = tmp_path / "odd.txt"
    odd.write_text("a\tb\tc\td\te\tf\n")
    with pytest.raises(ValueError, match="Could not tell which classifier"):
        detect_format(odd)


def test_selftest_passes(capsys):
    assert main(["selftest"]) == 0
    message = capsys.readouterr().err
    assert "is working" in message
    assert "5 of the 5 insertions were found" in message


def test_library_api(sim):
    import hgt_site

    sites = hgt_site.call_sites(f"{sim}.sam", [f"{sim}.kraken2.txt"], host_taxid=HOST_TAXID)
    good = hgt_site.passed(sites)
    assert len(good) == 5
    assert len(sites) > len(good)  # rejected candidates are returned too
    # Thresholds can be passed through by name, the same ones the CLI exposes.
    strict = hgt_site.passed(
        hgt_site.call_sites(f"{sim}.sam", [f"{sim}.kraken2.txt"], HOST_TAXID, min_support=1000)
    )
    assert strict == []
    pandas = pytest.importorskip("pandas")
    frame = hgt_site.to_dataframe(good)
    assert isinstance(frame, pandas.DataFrame)
    from hgt_site.call import COLUMNS

    assert list(frame.columns) == COLUMNS
    assert len(frame) == 5


def _call_to(prefix, out):
    assert main([
        "call", "--bam", f"{prefix}.sam", "--kraken", f"{prefix}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out),
    ]) == 0
    return out


def test_cohort_flags_sites_shared_by_every_sample(sim, tmp_path):
    # The same sample three times: every site is in all three, which is what a
    # protocol artefact looks like.
    table = _call_to(sim, tmp_path / "one.tsv")
    copies = []
    for name in "ABC":
        copy = tmp_path / f"{name}.tsv"
        copy.write_text(table.read_text())
        copies.append(str(copy))
    out, matrix = tmp_path / "cohort.tsv", tmp_path / "matrix.tsv"
    assert main(["cohort", *copies, "-o", str(out), "--matrix", str(matrix)]) == 0

    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 5
    for row in rows:
        assert row["n_samples"] == "3"
        assert row["n_samples_passed"] == "3"
        assert "COMMON" in row["flags"] and "RECURRENT" in row["flags"]
        assert row["samples"] == "A;B;C"
    header, *grid = matrix.read_text().splitlines()
    assert header.split("\t") == ["site", "A", "B", "C"]
    assert all(line.split("\t")[1:] == ["1", "1", "1"] for line in grid)


def test_cohort_keeps_sites_apart_when_samples_differ(tmp_path, tmp_path_factory):
    # Different seeds put insertions in different places, so nothing recurs.
    tables = []
    for seed in (1, 2):
        prefix = tmp_path_factory.mktemp(f"s{seed}") / "sample"
        simulate(str(prefix), SimConfig(seed=seed))
        tables.append(str(_call_to(str(prefix), tmp_path / f"seed{seed}.tsv")))
    out = tmp_path / "cohort.tsv"
    assert main(["cohort", *tables, "-o", str(out), "--names", "first", "second"]) == 0
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 10
    assert all(row["n_samples"] == "1" for row in rows)
    assert all(row["flags"] == "." for row in rows)
    assert {row["samples"] for row in rows} == {"first", "second"}


def test_cohort_rejects_a_mismatched_name_list(sim, tmp_path, capsys):
    table = str(_call_to(sim, tmp_path / "one.tsv"))
    assert main(["cohort", table, table, "-o", str(tmp_path / "c.tsv"), "--names", "only-one"]) == 1
    assert "one name per file" in capsys.readouterr().err


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


def test_unmapped_fastq_holds_every_unmapped_read(sim, tmp_path):
    from hgt_site.prepare import write_unmapped_fastq

    fastq = tmp_path / "unmapped.fq"
    n = write_unmapped_fastq(f"{sim}.sam", fastq)
    lines = fastq.read_text().splitlines()
    assert len(lines) == 4 * n
    names = {line[1:] for line in lines[::4]}
    expected = set()
    with open(f"{sim}.sam") as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            cols = line.split("\t")
            if int(cols[1]) & 0x4:
                expected.add(cols[0])
    assert names == expected
    for seq, plus, qual in zip(lines[1::4], lines[2::4], lines[3::4]):
        assert plus == "+" and len(qual) == len(seq) > 0


def test_unmapped_fastq_matches_between_sam_and_bam(sim, tmp_path):
    from hgt_site.prepare import write_unmapped_fastq

    pysam = pytest.importorskip("pysam")
    bam = tmp_path / "sample.bam"
    pysam.sort("-o", str(bam), f"{sim}.sam")
    from_sam, from_bam = tmp_path / "a.fq", tmp_path / "b.fq"
    write_unmapped_fastq(f"{sim}.sam", from_sam)
    write_unmapped_fastq(str(bam), from_bam)

    def records(path):
        lines = path.read_text().splitlines()
        return {(name[1:], seq) for name, seq in zip(lines[::4], lines[1::4])}

    assert records(from_sam) == records(from_bam)


def _fake_kraken2(directory, kraken_output):
    """A stand-in for kraken2 that copies a known per-read file to --output."""
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / "kraken2"
    script.write_text(
        "#!/bin/sh\n"
        "while [ $# -gt 0 ]; do\n"
        '  if [ "$1" = "--output" ]; then shift; out=$1; fi\n'
        "  shift\n"
        "done\n"
        f'cp "{kraken_output}" "$out"\n'
        'echo "fake kraken2 ran" >&2\n'
    )
    script.chmod(0o755)
    return directory


def test_run_command_classifies_and_calls(sim, tmp_path, monkeypatch):
    bin_dir = _fake_kraken2(tmp_path / "bin", f"{sim}.kraken2.txt")
    monkeypatch.setenv("PATH", str(bin_dir), prepend=":")
    database = tmp_path / "db"
    database.mkdir()
    out, work = tmp_path / "s.tsv", tmp_path / "work"
    assert main([
        "run", "--bam", f"{sim}.sam", "--kraken-db", str(database), "--host-taxid", str(HOST_TAXID),
        "-o", str(out), "--workdir", str(work),
    ]) == 0
    with open(out) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert sum(r["filter"] == "PASS" for r in rows) == 5
    assert (work / "unmapped.fq").exists() and (work / "kraken2.txt").exists()


def test_run_command_says_what_is_missing(sim, tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    database = tmp_path / "db"
    database.mkdir()
    code = main([
        "run", "--bam", f"{sim}.sam", "--kraken-db", str(database), "--host-taxid", str(HOST_TAXID),
        "-o", str(tmp_path / "s.tsv"),
    ])
    assert code == 1
    assert "kraken2 was not found" in capsys.readouterr().err


def _filters(path):
    with open(path) as handle:
        return {(r["chrom"], r["start"]): r["filter"] for r in csv.DictReader(handle, delimiter="\t")}


def test_excluded_region_stops_a_site_passing(sim, tmp_path):
    chrom, position = _truth_positions(sim)[0]
    bed = tmp_path / "blacklist.bed"
    bed.write_text(
        "track name=artefacts\n"
        f"{chrom}\t{position - 2000}\t{position + 2000}\n"
    )
    plain, filtered = tmp_path / "plain.tsv", tmp_path / "filtered.tsv"
    base = ["call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
            "--host-taxid", str(HOST_TAXID)]
    assert main(base + ["-o", str(plain)]) == 0
    assert main(base + ["-o", str(filtered), "--exclude-regions", str(bed)]) == 0

    before, after = _filters(plain), _filters(filtered)
    covered = [k for k, v in before.items() if k[0] == chrom and abs(int(k[1]) - position) < 2000]
    assert covered, "the simulated insertion should give a candidate site to exclude"
    for key in covered:
        assert "EXCLUDED" in after[key]
        assert after[key] != "PASS"
    # Sites outside the region are untouched.
    for key, value in before.items():
        if key not in covered:
            assert after[key] == value


def test_excluded_region_elsewhere_changes_nothing(sim, tmp_path):
    bed = tmp_path / "other.bed"
    bed.write_text("chrUnplaced\t0\t100000\n")
    plain, filtered = tmp_path / "a.tsv", tmp_path / "b.tsv"
    base = ["call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
            "--host-taxid", str(HOST_TAXID)]
    assert main(base + ["-o", str(plain)]) == 0
    assert main(base + ["-o", str(filtered), "--exclude-regions", str(bed)]) == 0
    assert _filters(plain) == _filters(filtered)


def test_regions_merge_and_respect_half_open_bounds(tmp_path):
    import gzip as _gzip

    from hgt_site.regions import load_regions, overlaps

    bed = tmp_path / "r.bed.gz"
    with _gzip.open(bed, "wt") as handle:
        handle.write("#comment\n\nchr1\t100\t200\nchr1\t150\t300\nchr1\t1000\t1100\nchr2\t0\t50\n")
    regions = load_regions(bed)
    assert regions["chr1"] == ((100, 1000), (300, 1100)), "overlapping intervals should merge"

    assert overlaps(regions, "chr1", 250, 260)
    assert overlaps(regions, "chr1", 50, 101), "a site ending one base inside the region overlaps"
    assert not overlaps(regions, "chr1", 50, 100), "BED end is exclusive, so touching is not overlapping"
    assert not overlaps(regions, "chr1", 300, 400), "BED start begins at 300 only for the next region"
    assert not overlaps(regions, "chr1", 400, 900), "between two regions"
    assert overlaps(regions, "chr2", 10, 20)
    assert not overlaps(regions, "chr3", 10, 20), "a contig with no regions"


def test_malformed_exclude_region_says_which_line(tmp_path):
    from hgt_site.regions import load_regions

    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t100\t200\nchr1\t300\n")
    with pytest.raises(ValueError, match="line 2"):
        load_regions(bed)

    bed.write_text("chr1\t500\t400\n")
    with pytest.raises(ValueError, match="not a valid interval"):
        load_regions(bed)


def test_donor_mode_finds_the_real_insertions(sim, tmp_path):
    """A donor FASTA should stand in for the classifier output.

    The site boundaries are not expected to match the classifier run exactly.
    Clustering is taxon-blind, so in the classifier run a contaminant anchor
    that happens to land beside a real site widens it, while in donor mode
    that read has no call at all. What must match is which insertions are
    found.
    """
    from hgt_site.simulate import DONOR_TAXID

    classified, from_donor = tmp_path / "classified.tsv", tmp_path / "donor.tsv"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(classified),
    ]) == 0
    assert main([
        "call", "--bam", f"{sim}.sam", "--donor", f"{sim}.donor.fa",
        "--host-taxid", str(HOST_TAXID), "--donor-taxid", str(DONOR_TAXID),
        "-o", str(from_donor),
    ]) == 0

    def passed(path):
        with open(path) as handle:
            return [
                (r["chrom"], int(r["start"]), int(r["end"]))
                for r in csv.DictReader(handle, delimiter="\t") if r["filter"] == "PASS"
            ]

    truth = _truth_positions(sim)
    for label, sites in (("classifier", passed(classified)), ("donor", passed(from_donor))):
        assert len(sites) == len(truth), f"{label} should report one site per insertion"
        for chrom, position in truth:
            assert any(
                c == chrom and start - 500 <= position <= end + 500 for c, start, end in sites
            ), f"{label} missed the insertion at {chrom}:{position}"


def test_donor_mode_reports_the_taxid_it_was_given(sim, tmp_path):
    out = tmp_path / "s.tsv"
    assert main([
        "call", "--bam", f"{sim}.sam", "--donor", f"{sim}.donor.fa",
        "--host-taxid", str(HOST_TAXID), "-o", str(out),
    ]) == 0
    with open(out) as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["filter"] == "PASS"]
    assert rows
    # The default, NCBI's taxon for a built sequence, when none was given.
    assert {r["top_taxid"] for r in rows} == {"32630"}


def test_donor_and_classifications_are_alternatives(sim, tmp_path, capsys):
    base = ["call", "--bam", f"{sim}.sam", "--host-taxid", str(HOST_TAXID),
            "-o", str(tmp_path / "s.tsv")]
    assert main(base) == 1, "one of the two must be given"
    assert "either --classifications" in capsys.readouterr().err
    both = base + ["--kraken", f"{sim}.kraken2.txt", "--donor", f"{sim}.donor.fa"]
    assert main(both) == 1, "giving both is ambiguous"
    assert "but not both" in capsys.readouterr().err


def test_donor_matching_is_strand_agnostic_and_needs_several_kmers():
    from hgt_site.donor import KMER, build_index, classify, revcomp

    donor = "".join(
        "ACGT"[(i * 7 + i // 3) % 4] for i in range(2000)
    )
    index = build_index({"d": donor})
    read = donor[500:650]
    calls = classify({"fwd": read, "rev": revcomp(read)}, index, taxid=7)
    assert calls == {"fwd": 7, "rev": 7}, "both orientations should match"

    # A read sharing less than the minimum number of k-mers is not called.
    shared = donor[:KMER + 2]
    unrelated = "".join("TTAGC"[(i * 3) % 5] for i in range(150 - len(shared)))
    assert classify({"edge": shared + unrelated}, index, taxid=7) == {}, (
        "a few k-mers of overlap should not be enough"
    )
    assert classify({"short": "ACGT"}, index, taxid=7) == {}, "a read shorter than k"


def test_donor_fasta_must_hold_sequence(tmp_path):
    from hgt_site.donor import read_fasta

    empty = tmp_path / "empty.fa"
    empty.write_text("")
    with pytest.raises(ValueError, match="no FASTA records"):
        read_fasta(empty)


def test_local_background_leaves_even_coverage_alone(sim, tmp_path):
    """Where coverage is even the option should change nothing but the column."""
    plain, local = tmp_path / "plain.tsv", tmp_path / "local.tsv"
    base = ["call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
            "--host-taxid", str(HOST_TAXID)]
    assert main(base + ["-o", str(plain)]) == 0
    assert main(base + ["-o", str(local), "--local-background"]) == 0

    def rows(path):
        with open(path) as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    before, after = rows(plain), rows(local)
    assert {r["depth_ratio"] for r in before} == {"1"}, "off by default"
    ratios = [float(r["depth_ratio"]) for r in after]
    assert any(r != 1.0 for r in ratios), "coverage should have been measured"
    assert all(0.5 < r < 2.0 for r in ratios), "even coverage should give ratios near 1"

    passed = lambda rs: {(r["chrom"], r["start"]) for r in rs if r["filter"] == "PASS"}
    assert passed(after) == passed(before)


def test_local_background_falls_back_on_a_filtered_alignment(sim, tmp_path, caplog):
    """Without the mapped pairs there is no coverage to measure."""
    filtered = tmp_path / "filtered.sam"
    with open(f"{sim}.sam") as source, open(filtered, "w") as out:
        for line in source:
            if line.startswith("@"):
                out.write(line)
                continue
            flag = int(line.split("\t")[1])
            if flag & 0x4 or flag & 0x8:
                out.write(line)

    out_tsv = tmp_path / "s.tsv"
    with caplog.at_level("WARNING"):
        assert main([
            "call", "--bam", str(filtered), "--kraken", f"{sim}.kraken2.txt",
            "--host-taxid", str(HOST_TAXID), "-o", str(out_tsv), "--local-background",
        ]) == 0
    assert "appears to have been filtered" in caplog.text
    with open(out_tsv) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows and {r["depth_ratio"] for r in rows} == {"1"}, "every site falls back"


def test_coverage_rate_and_guards():
    from hgt_site.alignments import Coverage
    from hgt_site.background import MAX_RATIO, MIN_LOCAL_READS, apply_ratios, unusable_reason
    from hgt_site.sites import Site

    coverage = Coverage(bin_size=1000)
    # chr1 bin 0 is busy, bin 1 is quiet, chr1 bin 5 is extreme.
    for _ in range(MIN_LOCAL_READS * 2):
        coverage.add(0, "chr1", 10)
    for _ in range(MIN_LOCAL_READS // 4):
        coverage.add(0, "chr1", 1010)
    for _ in range(MIN_LOCAL_READS * 100):
        coverage.add(0, "chr1", 5010)

    rate, reads = coverage.rate("chr1", 0, 1)
    assert reads == MIN_LOCAL_READS * 2
    assert rate == reads / 1000, "a narrow range still uses the whole bin"
    assert coverage.rate("chr1", 0, 1001)[1] == MIN_LOCAL_READS * 2 + MIN_LOCAL_READS // 4

    genome = 100_000
    busy, quiet, extreme = (Site("chr1", 0, 1), Site("chr1", 1000, 1001), Site("chr1", 5000, 5001))
    scaled = apply_ratios([busy, quiet, extreme], coverage, genome, window=0)
    assert busy.depth_ratio > 1, "a region above the mean should raise the expectation"
    assert quiet.depth_ratio == 1.0, "too few reads nearby to trust the estimate"
    assert extreme.depth_ratio == 1.0, f"a ratio above {MAX_RATIO} should be ignored"
    assert scaled == 1

    assert unusable_reason(None, genome) is not None
    assert unusable_reason(Coverage(bin_size=1000), genome) is not None, "no reads at all"


def test_unmapped_reads_do_not_count_as_coverage():
    from hgt_site.alignments import Coverage

    coverage = Coverage(bin_size=100)
    coverage.add(0x4, "chr1", 50)
    coverage.add(0, "*", -1)
    assert coverage.mapped_reads == 0 and not coverage.bins
    coverage.add(0x8, "chr1", 50)
    assert coverage.mapped_reads == 1, "a mapped read whose mate is unmapped is still coverage"
    assert coverage.fully_mapped_reads == 0, "but it is not a fully mapped pair"


def test_local_background_never_lowers_the_bar():
    """A quiet region is tested against the genome-wide rate, not a lower one.

    The ratio carries a few per cent of noise, so letting it fall below 1
    would lower the expected count and pass marginal candidates. On an even
    coverage benchmark that cost a false site.
    """
    from hgt_site.alignments import Coverage
    from hgt_site.background import MIN_LOCAL_READS, apply_ratios
    from hgt_site.sites import Site

    coverage = Coverage(bin_size=1000)
    for _ in range(MIN_LOCAL_READS * 4):
        coverage.add(0, "chr1", 10)
    for _ in range(MIN_LOCAL_READS):
        coverage.add(0, "chr1", 1010)

    # Genome mean sits between the two bins' rates.
    genome = 2000
    busy, quiet = Site("chr1", 0, 1), Site("chr1", 1000, 1001)
    apply_ratios([busy, quiet], coverage, genome, window=0, contig_lengths={"chr1": 2000})
    assert busy.depth_ratio > 1, "a crowded region should raise the expectation"
    assert quiet.depth_ratio == 1.0, "a quiet region must not lower it"


def test_composition_reports_how_foreign_the_reads_look(sim, tmp_path):
    """A real donor differs from the host in GC and in k-mer usage."""
    out = tmp_path / "s.tsv"
    assert main([
        "call", "--bam", f"{sim}.sam", "--kraken", f"{sim}.kraken2.txt",
        "--host-taxid", str(HOST_TAXID), "-o", str(out),
    ]) == 0
    with open(out) as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t") if r["filter"] == "PASS"]
    assert rows
    for row in rows:
        # The simulated donor is built at a higher GC than the host, so the
        # difference is positive and the k-mer profile sits well outside what
        # sampling noise would give.
        assert float(row["gc_difference"]) > 0.05, row["gc_difference"]
        assert float(row["kmer_ratio"]) > 1.2, row["kmer_ratio"]


def test_host_reads_called_foreign_look_like_the_host(sim, tmp_path):
    """The case the composition columns are actually for.

    Composition cannot separate a chance chimera from a real insertion, since
    both carry genuinely foreign sequence. What it can catch is a classifier
    calling host reads foreign, which leaves a site whose reads look exactly
    like the host.
    """
    import gzip
    import random

    # Build a classification file that calls the unmapped mates host-derived
    # sequence instead of the donor: same read names, but the sequence they
    # carry is host, so the composition should match the host.
    rng = random.Random(11)
    host = {}
    for line in open(f"{sim}.host.fa"):
        if line.startswith(">"):
            name = line[1:].split()[0]
            host[name] = []
        else:
            host[name].append(line.strip())
    host = {k: "".join(v) for k, v in host.items()}
    chrom, seq = next(iter(host.items()))

    sam_out = tmp_path / "hostlike.sam"
    calls_out = tmp_path / "hostlike.kraken2.txt"
    kept = []
    with open(f"{sim}.sam") as source, open(sam_out, "w") as sam:
        for line in source:
            if line.startswith("@"):
                sam.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            flag = int(f[1])
            if flag & 0x4:
                # Replace the unmapped mate's sequence with host sequence.
                start = rng.randint(0, len(seq) - len(f[9]))
                f[9] = seq[start : start + len(f[9])]
                kept.append(f[0])
            sam.write("\t".join(f) + "\n")
    with open(calls_out, "w") as handle:
        for name in kept:
            handle.write(f"C\t{name}\t562\t100\t562:50\n")

    out = tmp_path / "hostlike.tsv"
    assert main([
        "call", "--bam", str(sam_out), "--kraken", str(calls_out),
        "--host-taxid", str(HOST_TAXID), "-o", str(out), "--min-complexity", "0",
    ]) == 0
    import statistics

    with open(out) as handle:
        rows = [
            r for r in csv.DictReader(handle, delimiter="\t")
            if r["kmer_ratio"] and int(r["n_pairs"]) >= 5
        ]
    assert rows, "the run should produce candidate sites with several reads each"
    ratios = [float(r["kmer_ratio"]) for r in rows]
    differences = [abs(float(r["gc_difference"])) for r in rows]
    # A ratio near 1 is what sampling noise alone gives, which is the point:
    # these reads are host sequence wearing a foreign label. Individual sites
    # with a handful of reads are noisy, so this is about the typical site.
    assert statistics.median(ratios) < 1.5, f"host-derived reads look foreign: {ratios}"
    assert statistics.median(differences) < 0.05, f"and their GC is off: {differences}"
