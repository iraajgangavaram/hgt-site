"""Print the total Bioconda download count for hgt-site.

The Anaconda API has reported this in more than one shape over the years, so
rather than depend on one key this tries each in turn and says what it found
if none of them are there. A silent wrong number in the README is worse than a
failed workflow run.

    python scripts/fetch_conda_downloads.py

Prints a single integer on stdout. Anything else goes to stderr.
"""

import json
import sys
import urllib.error
import urllib.request

URL = "https://api.anaconda.org/package/bioconda/hgt-site"
TIMEOUT = 30


def total(payload):
    """The download count, from whichever shape this response uses."""
    for key in ("total_ndownloads", "ndownloads"):
        value = payload.get(key)
        if isinstance(value, int):
            return value
    files = payload.get("files")
    if isinstance(files, list) and files:
        counts = [f.get("ndownloads", 0) for f in files if isinstance(f, dict)]
        if any(isinstance(c, int) for c in counts):
            return sum(c for c in counts if isinstance(c, int))
    return None


def main():
    try:
        with urllib.request.urlopen(URL, timeout=TIMEOUT) as response:
            payload = json.load(response)
    except (urllib.error.URLError, json.JSONDecodeError, TimeoutError) as error:
        print(f"could not read {URL}: {error}", file=sys.stderr)
        return 1

    count = total(payload)
    if count is None:
        print(
            "no download count in the response. Top-level keys were: "
            + ", ".join(sorted(payload)),
            file=sys.stderr,
        )
        return 1
    print(count)
    return 0


if __name__ == "__main__":
    sys.exit(main())
