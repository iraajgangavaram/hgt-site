"""Produce the Kraken2 input that calling needs, straight from a host alignment.

hgt-site itself only needs a host alignment and Kraken2's per-read calls for
the reads that did not map. Making those calls takes two steps that are easy
to get wrong: pulling out exactly the unmapped reads, and keeping the read
names that tie them back to the alignment. This module does both, so a user
needs nothing but the alignment and a Kraken2 database.
"""

import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("hgt_site")

UNMAPPED, SECONDARY, QCFAIL, DUPLICATE, SUPPLEMENTARY = 0x4, 0x100, 0x200, 0x400, 0x800
SKIP = SECONDARY | SUPPLEMENTARY | QCFAIL | DUPLICATE


class ToolMissing(Exception):
    """A program hgt-site needs is not on PATH."""


def _fastq_record(name, seq, qual):
    if not qual or len(qual) != len(seq):
        qual = "I" * len(seq)
    return f"@{name}\n{seq}\n+\n{qual}\n"


def write_unmapped_fastq(alignments, path, reference=None, threads=1):
    """Write every unmapped read of a host alignment to a FASTQ file.

    Read names are kept exactly as the alignment holds them, so Kraken2's
    output lines up with the alignment. Both mates of a pair can appear, and
    a mate suffix is not added, because hgt-site matches on the base name.
    Returns the number of reads written.
    """
    path = Path(path)
    written = 0
    if str(alignments).endswith(".sam"):
        with open(alignments) as handle, open(path, "w") as out:
            for line in handle:
                if line.startswith("@"):
                    continue
                cols = line.rstrip("\n").split("\t")
                flag = int(cols[1])
                if not flag & UNMAPPED or flag & SKIP:
                    continue
                out.write(_fastq_record(cols[0], cols[9], cols[10]))
                written += 1
        return written

    import pysam

    mode = "rc" if str(alignments).endswith(".cram") else "rb"
    bam = pysam.AlignmentFile(str(alignments), mode, reference_filename=reference, threads=threads)
    with open(path, "w") as out:
        for read in bam.fetch(until_eof=True):
            if not read.flag & UNMAPPED or read.flag & SKIP:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            out.write(_fastq_record(read.query_name, seq, read.qual))
            written += 1
    bam.close()
    return written


def run_kraken2(fastq, database, output, report=None, threads=1, extra=()):
    """Classify the reads in a FASTQ file with Kraken2.

    Returns the path of the per-read output. Raises ToolMissing if Kraken2 is
    not installed, with the command that installs it.
    """
    kraken2 = shutil.which("kraken2")
    if kraken2 is None:
        raise ToolMissing(
            "kraken2 was not found on PATH. Install it with "
            "'conda install -c conda-forge -c bioconda kraken2', or classify the unmapped "
            "reads yourself and pass the result to hgt-site call --kraken."
        )
    database = Path(database)
    if not database.is_dir():
        raise ToolMissing(f"The Kraken2 database {database} is not a directory.")
    command = [
        kraken2, "--db", str(database), "--threads", str(threads),
        "--output", str(output), *extra, str(fastq),
    ]
    if report:
        command[5:5] = ["--report", str(report)]
    log.info("Running Kraken2 on the unmapped reads.")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise ToolMissing(f"Kraken2 failed:\n{result.stderr.strip()}")
    for line in result.stderr.strip().splitlines():
        if line.strip():
            log.info("Kraken2: %s", line.strip())
    return Path(output)


def prepare(alignments, database, workdir, reference=None, threads=1, kraken_args=()):
    """Write the unmapped reads and classify them. Returns the Kraken2 output path."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    fastq = workdir / "unmapped.fq"
    n = write_unmapped_fastq(alignments, fastq, reference, threads)
    log.info("Wrote %d unmapped reads to %s.", n, fastq)
    if n == 0:
        raise ToolMissing(
            "The alignment holds no unmapped reads, so there is nothing to classify. "
            "Check that the file is the alignment before any filtering of unmapped reads."
        )
    return run_kraken2(
        fastq, database, workdir / "kraken2.txt", workdir / "kraken2.report", threads, kraken_args
    )
