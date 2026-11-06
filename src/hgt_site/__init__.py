"""Call candidate horizontal gene transfer sites from host alignments and read classifications.

Everything the command line does is available here, so hgt-site can be used
from a notebook or inside another program without shelling out:

    >>> from hgt_site import call_sites, passed, to_dataframe
    >>> sites = call_sites("sample.bam", ["classified.txt"], host_taxid=9606)
    >>> to_dataframe(passed(sites)).head()

The pieces are also importable on their own: `load_classifications` for reading
a classifier's per-read output, `simulate` for building a dataset with known
insertion sites, and the `write_*` functions for the output formats.
"""

__version__ = "0.4.0"

__all__ = [
    "__version__",
    "DEFAULTS",
    "Site",
    "call_sites",
    "load_classifications",
    "passed",
    "simulate",
    "site_filter",
    "to_dataframe",
    "write_bed",
    "write_multiqc",
    "write_table",
    "write_vcf",
]


def __getattr__(name):
    # Imported on demand so that `import hgt_site` stays fast and does not
    # pull in pysam before it is needed.
    if name in ("DEFAULTS", "call_sites", "write_bed", "write_multiqc", "write_table", "write_vcf"):
        from . import call

        return getattr(call, name)
    if name in ("Site", "site_filter"):
        from . import sites

        return getattr(sites, name)
    if name == "load_classifications":
        from .classifiers import load_classifications

        return load_classifications
    if name == "simulate":
        from .simulate import simulate

        return simulate
    if name == "passed":
        return passed
    if name == "to_dataframe":
        return to_dataframe
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def passed(sites):
    """The sites that passed every filter, in the order they were called."""
    from .sites import site_filter

    return [site for site in sites if site_filter(site) == "PASS"]


def to_dataframe(sites):
    """Put sites into a pandas DataFrame with the same columns as the TSV.

    Needs pandas, which hgt-site does not install for you.
    """
    import pandas

    from .call import COLUMNS, _row

    return pandas.DataFrame([_row(site) for site in sites], columns=COLUMNS)
