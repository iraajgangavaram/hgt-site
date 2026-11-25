import io
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fetch_conda_downloads  # noqa: E402


def test_download_count_prefers_total():
    assert fetch_conda_downloads.total({"total_ndownloads": 17, "ndownloads": 11}) == 17


def test_download_count_accepts_alternate_top_level_key():
    assert fetch_conda_downloads.total({"ndownloads": 11}) == 11


def test_download_count_sums_valid_file_counts():
    payload = {"files": [{"ndownloads": 4}, {"ndownloads": 6}, {"ndownloads": "bad"}, {}]}

    assert fetch_conda_downloads.total(payload) == 10


def test_download_count_ignores_invalid_and_negative_counts():
    assert fetch_conda_downloads.total({"total_ndownloads": True, "files": [{"ndownloads": -1}]}) is None


def test_download_count_rejects_non_object_response():
    assert fetch_conda_downloads.total([1, 2, 3]) is None


def test_main_reports_malformed_response_without_traceback(monkeypatch, capsys):
    response = io.BytesIO(json.dumps([1, 2, 3]).encode())
    monkeypatch.setattr(fetch_conda_downloads.urllib.request, "urlopen", lambda *args, **kwargs: response)

    assert fetch_conda_downloads.main() == 1
    captured = capsys.readouterr()
    assert "response JSON type was list" in captured.err
    assert "Traceback" not in captured.err


def test_main_prints_download_count(monkeypatch, capsys):
    response = io.BytesIO(json.dumps({"total_ndownloads": 17}).encode())
    monkeypatch.setattr(fetch_conda_downloads.urllib.request, "urlopen", lambda *args, **kwargs: response)

    assert fetch_conda_downloads.main() == 0
    assert capsys.readouterr().out == "17\n"
