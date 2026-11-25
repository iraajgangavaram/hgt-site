import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmark"))
from run_benchmark import build_sample, run_pipeline  # noqa: E402


def test_run_pipeline_passes_paths_literally(tmp_path):
    value = "path with spaces; $(not-a-command)"
    output = tmp_path / "pipeline output.txt"
    producer = [sys.executable, "-c", "import sys; print(sys.argv[1])", value]
    consumer = [sys.executable, "-c", "import shutil, sys; shutil.copyfileobj(sys.stdin, sys.stdout)"]

    with output.open("w") as handle:
        run_pipeline(producer, consumer, stdout=handle)

    assert output.read_text() == value + "\n"


def test_run_pipeline_reports_producer_failure():
    producer = [sys.executable, "-c", "raise SystemExit(7)"]
    consumer = [sys.executable, "-c", "pass"]

    with pytest.raises(subprocess.CalledProcessError) as error:
        run_pipeline(producer, consumer)

    assert error.value.returncode == 7
    assert error.value.cmd == producer


def test_build_sample_rejects_more_insertions_than_host_can_hold():
    args = SimpleNamespace(n_insertions=2, insert_max=100, insert_min=50)

    with pytest.raises(SystemExit, match="host can hold at most 1 insertions"):
        build_sample({"host": "A" * 20_000}, {"donor": "C" * 100}, args, None)


def test_build_sample_allows_zero_insertions_on_a_short_host():
    args = SimpleNamespace(n_insertions=0, insert_max=100, insert_min=50)

    _, sample, truth = build_sample({"host": "A" * 100}, {"donor": "C" * 100}, args, None)

    assert sample == "A" * 100
    assert truth == []
