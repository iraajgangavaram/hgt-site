import pytest

from hgt_site.cohort import merge


@pytest.mark.parametrize(
    ("total", "present", "fraction", "is_common"),
    [
        (3, 2, 0.8, False),
        (4, 3, 0.8, False),
        (5, 3, 0.5, True),
        (5, 2, 0.5, False),
    ],
)
def test_common_fraction_requires_at_least_the_requested_share(total, present, fraction, is_common):
    samples = {
        f"sample-{index}": (
            [{"chrom": "chr1", "start": "100", "end": "200", "filter": "PASS"}]
            if index < present
            else []
        )
        for index in range(total)
    }

    rows = merge(samples, common_fraction=fraction)

    assert len(rows) == 1
    assert ("COMMON" in rows[0]["flags"]) is is_common
