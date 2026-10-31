import argparse
import sys

from . import __version__
from .call import DEFAULTS, call_sites, write_table
from .sites import site_filter


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hgt-site",
        description="Call candidate horizontal gene transfer sites from host alignments and Kraken2 read calls.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    call = sub.add_parser("call", help="Call candidate sites")
    call.add_argument("--bam", required=True, help="Paired-end reads aligned to the host (BAM, CRAM or SAM)")
    call.add_argument("--kraken", required=True, help="Kraken2 per-read output for the unmapped reads")
    call.add_argument("--host-taxid", required=True, type=int, help="NCBI taxid of the host, e.g. 9606")
    call.add_argument("--nodes", help="NCBI nodes.dmp, so the whole host lineage is excluded")
    call.add_argument("--exclude-taxid", type=int, action="append", default=[], help="Extra taxid to ignore (repeatable)")
    call.add_argument("--reference", help="Reference FASTA, needed for CRAM input")
    call.add_argument("-o", "--output", required=True, help="Output TSV")
    for key, value in DEFAULTS.items():
        call.add_argument(f"--{key.replace('_', '-')}", type=type(value), default=None, help=f"default {value}")
    sim = sub.add_parser("simulate", help="Write a small simulated dataset with known insertion sites")
    sim.add_argument("-o", "--output-prefix", required=True)
    sim.add_argument("--seed", type=int, default=1)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "call":
        opts = {k: getattr(args, k) for k in DEFAULTS}
        sites = call_sites(
            args.bam, args.kraken, args.host_taxid, args.nodes, args.exclude_taxid, args.reference, **opts
        )
        write_table(sites, args.output)
        passed = sum(1 for s in sites if site_filter(s) == "PASS")
        print(f"{len(sites)} candidate sites, {passed} passed filters -> {args.output}", file=sys.stderr)
    elif args.command == "simulate":
        from .simulate import SimConfig, simulate

        truth = simulate(args.output_prefix, SimConfig(seed=args.seed))
        print(f"Wrote {args.output_prefix}.sam, .kraken2.txt and .truth.bed ({len(truth)} insertions)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
