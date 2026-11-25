"""Input validation for cohort tables."""

import pytest

from hgt_site.cohort import read_sites


def test_read_sites_requires_cohort_columns(tmp_path):
    table = tmp_path / "sample.tsv"
    table.write_text("chrom\tstart\tend\tn_pairs\n")

    with pytest.raises(ValueError, match=r"missing column\(s\).*filter"):
        read_sites(table)


@pytest.mark.parametrize(
    ("header", "row", "message"),
    [
        ("chrom\tstart\tend\tn_pairs\tfilter", "\t10\t20\t3\tPASS", "chromosome"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\tx\t20\t3\tPASS", "start must be an integer"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\t10\tnope\t3\tPASS", "end must be an integer"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\t10\t20\tnope\tPASS", "n_pairs must be an integer"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\t-1\t20\t3\tPASS", "start must be non-negative"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\t21\t20\t3\tPASS", "end must be greater than or equal to start"),
        ("chrom\tstart\tend\tn_pairs\tfilter", "chr1\t10\t20\t-1\tPASS", "n_pairs must be non-negative"),
    ],
)
def test_read_sites_reports_invalid_selected_rows_with_line_context(
    tmp_path, header, row, message
):
    table = tmp_path / "sample.tsv"
    table.write_text(f"{header}\n{row}\n")

    with pytest.raises(ValueError, match=message) as error:
        read_sites(table)

    assert f"{table}:2:" in str(error.value)


def test_read_sites_skips_invalid_failed_rows_unless_requested(tmp_path):
    table = tmp_path / "sample.tsv"
    table.write_text(
        "chrom\tstart\tend\tn_pairs\tfilter\n"
        "chr1\tbad\t20\t3\tLOW_SUPPORT\n"
        "chr2\t10\t20\t4\tPASS\n"
    )

    rows = read_sites(table)
    assert len(rows) == 1
    assert (rows[0]["chrom"], rows[0]["start"], rows[0]["end"], rows[0]["n_pairs"]) == (
        "chr2", 10, 20, 4
    )
    with pytest.raises(ValueError, match=r"sample\.tsv:2: start must be an integer"):
        read_sites(table, include_failed=True)


def test_read_sites_rejects_headerless_input(tmp_path):
    table = tmp_path / "empty.tsv"
    table.write_text("")

    with pytest.raises(ValueError, match="no header row"):
        read_sites(table)
