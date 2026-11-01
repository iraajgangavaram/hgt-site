import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .call import DEFAULTS, call_sites, write_bed, write_multiqc, write_table, write_vcf
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


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hgt-site",
        description="Call candidate horizontal gene transfer sites from host alignments and Kraken2 read calls.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser("call", help="Call candidate sites")
    call.add_argument("--bam", required=True, help="Paired-end reads aligned to the host (BAM, CRAM or SAM)")
    call.add_argument("--kraken", required=True, nargs="+", help="Kraken2 per-read output for the unmapped reads (one or more files)")
    call.add_argument("--host-taxid", required=True, type=int, help="NCBI taxid of the host, e.g. 9606")
    call.add_argument("--nodes", help="NCBI nodes.dmp, so the whole host lineage is excluded")
    call.add_argument("--exclude-taxid", type=int, action="append", default=[], help="Extra taxid to ignore (repeatable)")
    call.add_argument("--reference", help="Reference FASTA, needed for CRAM input")
    call.add_argument("--threads", type=int, default=1, help="Threads for reading BAM or CRAM files (default 1)")
    call.add_argument(
        "--no-split-reads", action="store_true",
        help="Do not look for split reads at the junctions (they are only used with indexed BAM or CRAM, or SAM)",
    )
    call.add_argument("-o", "--output", required=True, help="Output TSV")
    call.add_argument("--bed", help="Also write the sites that passed every filter as a BED file")
    call.add_argument("--vcf", help="Also write every candidate site as a VCF file, with PASS or the failing flags in FILTER")
    call.add_argument("--multiqc", help="Also write a summary table for MultiQC (name it <something>_mqc.json)")
    call.add_argument("--sample-name", help="Sample name for the MultiQC table (default: the alignment file name)")
    limits = call.add_argument_group("thresholds")
    for key, value in DEFAULTS.items():
        limits.add_argument(
            f"--{key.replace('_', '-')}", type=type(value), default=None,
            help=f"{THRESHOLD_HELP[key]} (default {value})",
        )
    sim = sub.add_parser("simulate", help="Write a small simulated dataset with known insertion sites")
    sim.add_argument("-o", "--output-prefix", required=True)
    sim.add_argument("--seed", type=int, default=1)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
    if args.command == "call":
        opts = {k: getattr(args, k) for k in DEFAULTS}
        summary = {}
        sites = call_sites(
            args.bam, args.kraken, args.host_taxid, args.nodes, args.exclude_taxid, args.reference,
            threads=args.threads, split_reads=not args.no_split_reads, summary=summary, **opts
        )
        write_table(sites, args.output)
        if args.bed:
            write_bed(sites, args.bed)
        if args.vcf:
            write_vcf(sites, args.vcf, __version__)
        if args.multiqc:
            sample = args.sample_name or Path(args.bam).name.split(".")[0]
            write_multiqc(sites, args.multiqc, sample, summary)
        passed = sum(1 for s in sites if site_filter(s) == "PASS")
        print(f"{len(sites)} candidate sites, {passed} passed filters -> {args.output}", file=sys.stderr)
    elif args.command == "simulate":
        from .simulate import SimConfig, simulate

        truth = simulate(args.output_prefix, SimConfig(seed=args.seed))
        print(f"Wrote {args.output_prefix}.sam, .kraken2.txt, .truth.bed and .host.fa ({len(truth)} insertions)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
