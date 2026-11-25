import pytest

from hgt_site.control import load_control


def test_control_rejects_missing_required_columns(tmp_path):
    table = tmp_path / "control.tsv"
    table.write_text("chrom\tstart\tend\nchr1\t10\t20\n")

    with pytest.raises(ValueError, match="missing required column.*n_pairs"):
        load_control(table)


def test_control_rejects_invalid_pair_count_instead_of_dropping_site(tmp_path):
    table = tmp_path / "control.tsv"
    table.write_text("chrom\tstart\tend\tn_pairs\nchr1\t10\t20\tnot-a-count\n")

    with pytest.raises(ValueError, match=r"control.tsv line 2: n_pairs must be a whole number"):
        load_control(table)


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("chr1\tbad\t20\t8", "start must be a whole number"),
        ("chr1\t10\tbad\t8", "end must be a whole number"),
        ("chr1\t20\t10\t8", "20-10 is not a valid interval"),
        ("chr1\t10\t20\t-1", "n_pairs must not be negative"),
        ("\t10\t20\t8", "chrom must not be empty"),
    ],
)
def test_control_reports_invalid_row_with_line_context(tmp_path, row, message):
    table = tmp_path / "control.tsv"
    table.write_text(f"chrom\tstart\tend\tn_pairs\n{row}\n")

    with pytest.raises(ValueError, match=rf"control.tsv line 2: {message}"):
        load_control(table)
