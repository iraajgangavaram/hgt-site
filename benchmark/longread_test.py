"""Benchmark hgt-site on long reads.

A long read that crosses an integration carries both sides of the junction in
one molecule, so the aligner places the host part and clips the rest. There is
no mate to infer anything from and the junction is exact rather than an
interval, which should make this the easier case.

This builds a host with known insertions, simulates long reads from it with a
nanopore-like error rate and length distribution, aligns them with minimap2 and
calls the sites with `--long-reads`.

    python benchmark/longread_test.py --host host.fa --donor ecoli.fa \\
        --contaminants klebsiella.fa --outdir out

The error model is a flat per-base substitution and indel rate rather than the
context-dependent errors a real basecaller makes, and the length distribution
is a simple gamma. It is enough to test whether the junctions are found and
whether chance chimeras are rejected, and not enough to say anything about how
hgt-site behaves on a particular chemistry.
"""

import argparse
import csv
import math
import random
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_benchmark import concat, evaluate, read_fasta, revcomp, run, write_fasta  # noqa: E402

BASES = "ACGT"


def noisy(seq, rng, sub=0.04, ins=0.015, dele=0.015):
    """Add nanopore-like noise: a few per cent substitutions and short indels."""
    out = []
    for base in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(rng.choice(BASES))
        if rng.random() < sub:
            out.append(rng.choice([b for b in BASES if b != base]))
        else:
            out.append(base)
    return "".join(out)


def read_lengths(n, rng, mean=12_000, shape=2.0):
    """Gamma-distributed read lengths, which is the usual shape for nanopore."""
    scale = mean / shape
    return [max(1_000, int(rng.gammavariate(shape, scale))) for _ in range(n)]


def simulate(seq, n, rng, out, prefix, mean_length):
    """Write n long reads sampled from seq."""
    written = 0
    with open(out, "w") as handle:
        for i, length in enumerate(read_lengths(n, rng, mean_length)):
            if length >= len(seq):
                continue
            start = rng.randint(0, len(seq) - length)
            piece = seq[start : start + length]
            if rng.random() < 0.5:
                piece = revcomp(piece)
            piece = noisy(piece, rng)
            handle.write(f"@{prefix}_{i}\n{piece}\n+\n{'I' * len(piece)}\n")
            written += 1
    return written


def build_sample(host_seq, donor_seq, n_insertions, rng, insert_min, insert_max):
    """Insert donor fragments into the host at well separated positions."""
    positions = []
    while len(positions) < n_insertions:
        pos = rng.randint(50_000, len(host_seq) - 50_000)
        if all(abs(pos - p) > 100_000 for p in positions):
            positions.append(pos)
    positions.sort()

    pieces, truth, cursor = [], [], 0
    for pos in positions:
        length = rng.randint(insert_min, insert_max)
        start = rng.randint(0, len(donor_seq) - length)
        fragment = donor_seq[start : start + length]
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        pieces += [host_seq[cursor:pos], fragment]
        truth.append(pos)
        cursor = pos
    pieces.append(host_seq[cursor:])
    return "".join(pieces), truth


def chimeras(host_seq, others, n, rng, out, mean_length):
    """Long chance chimeras: one long host molecule joined to a foreign one.

    Rarer in a real long-read library than in a short-read one, because two
    long molecules have to be joined rather than two short ones, but they do
    happen and they are the artefact the test is meant to reject.
    """
    half = mean_length // 2
    # Short contigs, such as plasmids, cannot supply half a long molecule.
    names = [n for n, s in others.items() if len(s) > half + 1]
    if not names:
        raise SystemExit(
            f"no contaminant contig is longer than {half + 1} bp, so no long chimera can be built"
        )
    with open(out, "w") as handle:
        for i in range(n):
            hpos = rng.randint(0, len(host_seq) - half - 1)
            source = others[rng.choice(names)]
            cpos = rng.randint(0, len(source) - half - 1)
            piece = host_seq[hpos : hpos + half] + source[cpos : cpos + half]
            if rng.random() < 0.5:
                piece = revcomp(piece)
            piece = noisy(piece, rng)
            handle.write(f"@chimera_{i}\n{piece}\n+\n{'I' * len(piece)}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--host", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--contaminants", nargs="+", default=[])
    p.add_argument("--outdir", required=True)
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--host-taxid", type=int, default=9606)
    p.add_argument("--donor-taxid", type=int, default=562)
    p.add_argument("--n-insertions", type=int, default=20)
    p.add_argument("--insert-min", type=int, default=2_000)
    p.add_argument("--insert-max", type=int, default=8_000)
    p.add_argument("--coverage", type=float, default=5.0)
    p.add_argument("--mean-length", type=int, default=12_000)
    p.add_argument("--contaminant-reads", type=int, default=2_000)
    p.add_argument("--chimeras", type=int, default=500)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    host_seq = max(read_fasta(args.host).values(), key=len).upper()
    donor_seq = max(read_fasta(args.donor).values(), key=len).upper()
    if len(donor_seq) < args.insert_max:
        raise SystemExit(f"donor is {len(donor_seq)} bp, shorter than --insert-max")

    reference = out / "host.fa"
    write_fasta(reference, {"host": host_seq})
    sample_seq, truth = build_sample(
        host_seq, donor_seq, args.n_insertions, rng, args.insert_min, args.insert_max
    )
    write_fasta(out / "sample.fa", {"sample": sample_seq})
    with open(out / "truth.bed", "w") as handle:
        for pos in truth:
            handle.write(f"host\t{pos}\t{pos + 1}\n")

    n_reads = math.ceil(args.coverage * len(sample_seq) / args.mean_length)
    print(f"host {len(host_seq):,} bp, {len(truth)} insertions, {n_reads:,} long reads")
    simulate(sample_seq, n_reads, rng, out / "sample.fq", "sample", args.mean_length)

    parts = [out / "sample.fq"]
    others = {}
    for i, path in enumerate(args.contaminants):
        others.update(read_fasta(path))
        seq = max(read_fasta(path).values(), key=len).upper()
        simulate(seq, args.contaminant_reads, rng, out / f"c{i}.fq", f"contam{i}", args.mean_length)
        parts.append(out / f"c{i}.fq")
    if others and args.chimeras:
        chimeras(host_seq, others, args.chimeras, rng, out / "chim.fq", args.mean_length)
        parts.append(out / "chim.fq")
    concat(parts, out / "all.fq")

    with open(out / "aligned.sam", "w") as handle:
        subprocess.run(
            [args.minimap2, "-ax", "map-ont", "-t", str(args.threads),
             str(reference), str(out / "all.fq")],
            stdout=handle, stderr=subprocess.DEVNULL, check=True,
        )
    run(["samtools", "sort", "-o", str(out / "aligned.bam"), str(out / "aligned.sam")])
    run(["samtools", "index", str(out / "aligned.bam")])

    run(["hgt-site", "call", "--bam", str(out / "aligned.bam"), "--long-reads",
         "--donor", args.donor, "--donor-taxid", str(args.donor_taxid),
         "--host-taxid", str(args.host_taxid), "--threads", str(args.threads),
         "-o", str(out / "sites.tsv")])

    with open(out / "sites.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    calls = [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows if r["filter"] == "PASS"]
    stats = evaluate(calls, [("host", pos, 0) for pos in truth])

    # How far the reported junction sits from the real one. A basecaller's
    # errors move the clip boundary by a few bases, so asking for the exact
    # base understates how well the junction is placed.
    offsets = sorted(
        min(abs(int(r["start"]) - t) for t in truth)
        for r in rows if r["filter"] == "PASS"
    )
    exact = sum(1 for o in offsets if o == 0)
    within10 = sum(1 for o in offsets if o <= 10)
    print(f"\ncandidates      {len(rows)}")
    print(f"sites reported  {stats['calls']}")
    print(f"of those real   {stats['true_calls']}")
    print(f"precision       {stats['precision']:.3f}")
    print(f"recall          {stats['recall']:.3f} ({stats['recovered']} of {stats['truth']})")
    if offsets:
        print(f"junction exact   {exact} of {len(offsets)}")
        print(f"within 10 bp     {within10} of {len(offsets)}")
        print(f"median offset    {offsets[len(offsets) // 2]} bp")
    print("\nfilter column:")
    for flag, n in Counter(r["filter"] for r in rows).most_common(6):
        print(f"  {flag:<26} {n}")


if __name__ == "__main__":
    main()
