import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .call import DEFAULTS, call_sites, write_bed, write_multiqc, write_table, write_vcf
from .classifiers import FORMATS
from .cohort import COMMON_FRACTION
from .cohort import WINDOW as COHORT_WINDOW
from .sites import site_filter

THRESHOLD_HELP = {
    "min_mapq": "Minimum mapping quality for a read to count as an anchor",
    "max_gap": "Maximum gap in bp between neighbouring anchors in one site",
    "max_span": "Maximum width in bp of one site",
    "window": "Window in bp used for the enrichment test and the split-read search",
    "min_support": "Minimum number of supporting pairs for a site to pass",
    "min_taxon_fraction": "Minimum share of pairs from the main taxon before MIXED_TAXA is considered",
    "max_other_pvalue": "P-value below which pairs from other taxa count as an excess (MIXED_TAXA)",
    "min_mean_mapq": "Minimum mean mapping quality of the anchors (LOW_MAPQ)",
    "min_complexity": "Minimum mean 3-mer entropy of the unmapped mates (LOW_COMPLEXITY)",
    "max_qvalue": "Maximum q-value of the enrichment test (NOT_ENRICHED)",
}


def _add_common(parser, kraken_required):
    parser.add_argument("--bam", required=True, help="Paired-end reads aligned to the host (BAM, CRAM or SAM)")
    if kraken_required:
        parser.add_argument(
            "--classifications", "--kraken", dest="kraken", required=True, nargs="+",
            help="Per-read taxonomic calls for the unmapped reads, from Kraken2, KrakenUniq, Kaiju, "
                 "Centrifuge, DIAMOND or BLAST, or any two-column file of read name and taxid "
                 "(one or more files, gzipped or not)",
        )
    parser.add_argument("--host-taxid", required=True, type=int, help="NCBI taxid of the host, e.g. 9606")
    parser.add_argument("--nodes", help="NCBI nodes.dmp, so the whole host lineage is excluded")
    parser.add_argument(
        "--exclude-taxid", type=int, action="append", default=[], help="Extra taxid to ignore (repeatable)"
    )
    parser.add_argument(
        "--classifier", default="auto", choices=FORMATS,
        help="Which classifier wrote the per-read calls (default: work it out from the file)",
    )
    parser.add_argument("--reference", help="Reference FASTA, needed for CRAM input")
    parser.add_argument(
        "--annotation",
        help="GTF or GFF3 file of gene annotations, gzipped or not. Each site is reported with the gene "
             "it falls in, or the nearest gene within --max-gene-distance",
    )
    parser.add_argument(
        "--max-gene-distance", type=int, default=100_000, metavar="BP",
        help="How far from a gene a site can be and still have it reported as the nearest (default 100000)",
    )
    parser.add_argument("--threads", type=int, default=1, help="Threads for reading files and for Kraken2 (default 1)")
    parser.add_argument(
        "--no-split-reads", action="store_true",
        help="Do not look for split reads at the junctions (they are only used with indexed BAM or CRAM, or SAM)",
    )
    parser.add_argument("-o", "--output", required=True, help="Output TSV")
    parser.add_argument("--bed", help="Also write the sites that passed every filter as a BED file")
    parser.add_argument(
        "--vcf",
        help="Also write every candidate site as a VCF file, with PASS or the failing flags in FILTER. "
             "A name ending in .gz is block compressed and indexed with tabix",
    )
    parser.add_argument("--multiqc", help="Also write a summary table for MultiQC (name it <something>_mqc.json)")
    parser.add_argument("--sample-name", help="Sample name for the MultiQC table (default: the alignment file name)")
    limits = parser.add_argument_group("thresholds")
    for key, value in DEFAULTS.items():
        limits.add_argument(
            f"--{key.replace('_', '-')}", type=type(value), default=None,
            help=f"{THRESHOLD_HELP[key]} (default {value})",
        )


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hgt-site",
        description="Call candidate horizontal gene transfer sites from host alignments and Kraken2 read calls.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser(
        "run",
        help="Classify the unmapped reads with Kraken2 and call sites, in one step",
        description="Pull the unmapped reads out of a host alignment, classify them with Kraken2, and call "
                    "candidate sites. Needs kraken2 on PATH and a Kraken2 database.",
    )
    _add_common(run, kraken_required=False)
    run.add_argument("--kraken-db", required=True, help="Kraken2 database directory")
    run.add_argument(
        "--workdir",
        help="Where to keep the unmapped reads and the Kraken2 output (default: a temporary folder, deleted after)",
    )
    run.add_argument(
        "--kraken-arg", action="append", default=[], metavar="ARG",
        help="Extra argument passed on to Kraken2, such as --confidence=0.1 (repeatable)",
    )

    call = sub.add_parser("call", help="Call candidate sites from an alignment and existing Kraken2 output")
    _add_common(call, kraken_required=True)

    cohort = sub.add_parser(
        "cohort",
        help="Combine the sites called in several samples and flag the recurrent ones",
        description="Merge the per-sample tables hgt-site writes into one, so a position seen in "
                    "every sample, which is usually an artefact, can be told from one seen in a "
                    "few. Reads the output of calling, so samples called at different times can "
                    "be combined without rerunning anything.",
    )
    cohort.add_argument("sites", nargs="+", help="Per-sample hgt-site TSV files, one per sample")
    cohort.add_argument(
        "--names", nargs="+",
        help="Sample names, in the same order as the files (default: taken from the file names)",
    )
    cohort.add_argument("-o", "--output", required=True, help="Output TSV of merged positions")
    cohort.add_argument("--matrix", help="Also write a position by sample presence grid")
    cohort.add_argument(
        "--window", type=int, default=COHORT_WINDOW, metavar="BP",
        help=f"Positions this far apart or closer count as the same site (default {COHORT_WINDOW})",
    )
    cohort.add_argument(
        "--common-fraction", type=float, default=COMMON_FRACTION, metavar="F",
        help=f"Flag a site as COMMON when it appears in this share of the samples (default {COMMON_FRACTION})",
    )
    cohort.add_argument(
        "--include-failed", action="store_true",
        help="Use every candidate site, not only the ones that passed every filter",
    )

    sub.add_parser(
        "selftest",
        help="Check the installation by calling a dataset whose answers are known",
        description="Build a small dataset with known insertion sites, call it, and check the "
                    "sites that come back are the ones that were put in. Needs no data and no "
                    "database, and leaves nothing behind.",
    )

    sim = sub.add_parser("simulate", help="Write a small simulated dataset with known insertion sites")
    sim.add_argument("-o", "--output-prefix", required=True)
    sim.add_argument("--seed", type=int, default=1)
    return parser


def _write_outputs(args, sites, summary):
    write_table(sites, args.output)
    sample = args.sample_name or Path(args.bam).name.split(".")[0]
    if args.bed:
        write_bed(sites, args.bed)
    if args.vcf:
        write_vcf(
            sites, args.vcf, __version__,
            contigs=summary.get("contig_lengths"),
            command=args.command_line,
            settings=summary.get("settings"),
            sample=sample,
            reference=args.reference,
        )
    if args.multiqc:
        write_multiqc(sites, args.multiqc, sample, summary)
    passed = sum(1 for s in sites if site_filter(s) == "PASS")
    print(f"{len(sites)} candidate sites, {passed} passed filters -> {args.output}", file=sys.stderr)


def _call(args, kraken):
    opts = {k: getattr(args, k) for k in DEFAULTS}
    summary = {}
    sites = call_sites(
        args.bam, kraken, args.host_taxid, args.nodes, args.exclude_taxid, args.reference,
        threads=args.threads, split_reads=not args.no_split_reads, annotation=args.annotation,
        max_gene_distance=args.max_gene_distance, classifier=args.classifier, summary=summary, **opts
    )
    _write_outputs(args, sites, summary)


def _selftest():
    """Call a dataset with known answers and report whether they came back.

    This exists so that someone can tell a broken installation from a broken
    input before they spend time on their own data.
    """
    import tempfile

    from .evaluate import evaluate, read_truth
    from .simulate import SimConfig, simulate

    with tempfile.TemporaryDirectory(prefix="hgt-site-selftest-") as folder:
        prefix = str(Path(folder) / "selftest")
        truth = simulate(prefix, SimConfig(seed=1))
        print(f"Built a dataset with {len(truth)} known insertions.", file=sys.stderr)
        sites = call_sites(f"{prefix}.sam", [f"{prefix}.kraken2.txt"], 9606)
        passed = [s for s in sites if site_filter(s) == "PASS"]
        calls = [(s.chrom, s.start, s.end) for s in passed]
        score = evaluate(calls, read_truth(f"{prefix}.truth.bed"))

    print(
        f"Called {len(sites)} candidate sites, of which {len(passed)} passed every filter.",
        file=sys.stderr,
    )
    print(
        f"Of those, {score['true_calls']} are real insertions, and "
        f"{score['recovered']} of the {score['truth']} insertions were found.",
        file=sys.stderr,
    )
    if score["true_calls"] == len(passed) and score["recovered"] == score["truth"]:
        print(f"hgt-site {__version__} is working.", file=sys.stderr)
        return 0
    print(
        "This installation did not reproduce the expected result. Please report it at "
        "https://github.com/iraajgangavaram/hgt-site/issues with this output and your "
        "Python version.",
        file=sys.stderr,
    )
    return 1


def main(argv=None):
    given = sys.argv[1:] if argv is None else list(argv)
    args = build_parser().parse_args(given)
    # Recorded in the VCF header. Taken from what was parsed rather than from
    # sys.argv, so an embedded call records its own arguments and not the
    # arguments of whatever program called it.
    args.command_line = " ".join(["hgt-site", *given])
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
    if args.command == "call":
        _call(args, args.kraken)
    elif args.command == "run":
        import tempfile

        from .prepare import ToolMissing, prepare

        temporary = None
        if args.workdir:
            workdir = args.workdir
        else:
            temporary = tempfile.TemporaryDirectory(prefix="hgt-site-")
            workdir = temporary.name
        try:
            kraken = prepare(
                args.bam, args.kraken_db, workdir, args.reference, args.threads, args.kraken_arg
            )
            _call(args, [str(kraken)])
        except ToolMissing as problem:
            print(f"ERROR: {problem}", file=sys.stderr)
            return 1
        finally:
            if temporary is not None:
                temporary.cleanup()
    elif args.command == "cohort":
        from .cohort import build, write_cohort, write_matrix

        try:
            rows, names = build(
                args.sites, args.names, args.window, args.common_fraction, args.include_failed
            )
        except ValueError as problem:
            print(f"ERROR: {problem}", file=sys.stderr)
            return 1
        write_cohort(rows, args.output)
        if args.matrix:
            write_matrix(rows, names, args.matrix)
        common = sum(1 for r in rows if "COMMON" in r["flags"])
        print(
            f"{len(rows)} positions across {len(names)} samples, {common} of them in most samples "
            f"-> {args.output}",
            file=sys.stderr,
        )
    elif args.command == "selftest":
        return _selftest()
    elif args.command == "simulate":
        from .simulate import SimConfig, simulate

        truth = simulate(args.output_prefix, SimConfig(seed=args.seed))
        print(
            f"Wrote {args.output_prefix}.sam, .kraken2.txt, .truth.bed and .host.fa ({len(truth)} insertions)",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
