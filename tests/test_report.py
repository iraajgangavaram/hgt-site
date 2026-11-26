"""Tests for the standalone results report."""

from hgt_site.cli import main


HEADER = "chrom\tstart\tend\tn_pairs\ttop_taxid\tfilter\tflags\tqvalue\tgene\tleft_breakpoint\tleft_split_reads\tright_breakpoint\tright_split_reads\n"


def test_report_command_writes_searchable_self_contained_html(tmp_path, capsys):
    sites = tmp_path / "sample.hgt-site.tsv"
    sites.write_text(
        HEADER
        + "chr1\t100\t200\t5\t562\tPASS\tONE_SIDED\t0.001\tgeneA\t150\t2\t\t0\n"
        + "chr2\t300\t400\t2\t9606\tLOW_SUPPORT;NOT_ENRICHED\tLOW_SUPPORT;NOT_ENRICHED\t0.8\tgeneB\t\t0\t\t0\n"
    )
    report = tmp_path / "report.html"

    assert main(["report", str(sites), "-o", str(report)]) == 0
    output = report.read_text()

    assert "2 candidate sites" in capsys.readouterr().err
    assert "<strong>2</strong><span>candidate sites" in output
    assert "<strong>1</strong><span>passed filters" in output
    assert "<strong>1</strong><span>did not pass" in output
    assert "<strong>1</strong><span>with split-read junctions" in output
    assert "LOW_SUPPORT" in output and "NOT_ENRICHED" in output
    assert "geneA" in output and "geneB" in output
    assert "Estimated from read pairs" in output
    assert 'type="search"' in output
    assert "PAGE_SIZE = 50" in output
    assert 'window.addEventListener("beforeprint", () => render(true))' in output
    assert 'href="sample.hgt-site.tsv"' in output
    assert "<script src=" not in output
    assert "<link rel=\"stylesheet\"" not in output


def test_report_escapes_untrusted_table_values(tmp_path):
    sites = tmp_path / "unsafe.tsv"
    sites.write_text(HEADER + "chr1\t1\t2\t4\t562\tPASS\t.\t0.01\t</script><img src=x onerror=alert(1)>\t\t0\t\t0\n")
    report = tmp_path / "unsafe.html"

    assert main(["report", str(sites), "-o", str(report)]) == 0
    output = report.read_text()

    assert "</script><img" not in output
    assert "\\u003c/script>" in output


def test_report_accepts_a_header_only_table(tmp_path):
    sites = tmp_path / "empty.tsv"
    sites.write_text(HEADER)
    report = tmp_path / "empty.html"

    assert main(["report", str(sites), "-o", str(report)]) == 0
    output = report.read_text()

    assert "<strong>0</strong><span>candidate sites" in output
    assert "No candidate sites match this search." in output


def test_report_rejects_incomplete_input_without_writing_output(tmp_path, capsys):
    sites = tmp_path / "incomplete.tsv"
    sites.write_text("chrom\tstart\tend\n")
    report = tmp_path / "report.html"

    assert main(["report", str(sites), "-o", str(report)]) == 1
    error = capsys.readouterr().err
    assert error.startswith("ERROR:")
    assert "missing required column(s)" in error
    assert "Traceback" not in error
    assert not report.exists()


def test_report_cannot_overwrite_its_input(tmp_path, capsys):
    sites = tmp_path / "sample.tsv"
    sites.write_text(HEADER)

    assert main(["report", str(sites), "-o", str(sites)]) == 1
    assert "would overwrite input site table" in capsys.readouterr().err
    assert sites.read_text() == HEADER
