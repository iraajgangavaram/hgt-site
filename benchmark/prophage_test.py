"""Benchmark a bacterial host with an element inserted into a repeated gene.

The question this answers came from someone asking whether hgt-site can be used
on prophage, ICE and IME integrations. The method does not care what the host
is, so a bacterial host should work, but prophages and ICEs characteristically
integrate into tRNA genes, and a tRNA gene family is a set of near-identical
copies scattered around the genome. Host reads landing on one copy cannot be
placed confidently, so the anchors that hold a site up may fall below the
mapping quality cut and the site is reported as LOW_MAPQ rather than passing.

That is a prediction about the method rather than an observation, so this
measures it. Two conditions over the same host and the same element:

    adjacent  each insertion sits immediately after a copy of a repeated
              300 bp sequence, as an integration into a tRNA gene does
    unique    each insertion sits in sequence that occurs once, which is the
              condition every other benchmark on this tool has used

Both run against a host that carries the same repeat family, so the only thing
that changes between them is whether the insertion is inside it.

    python benchmark/prophage_test.py --host ecoli.fa --element pseudomonas.fa \\
        --contaminants klebsiella.fa salmonella.fa \\
        --kraken-db k2db --nodes k2db/taxonomy/nodes.dmp --outdir out

The element stands in for a prophage by being a slice of a genome from another
genus: distinct base composition, and its own taxid in the classifier, which is
what a real prophage gives you. It is not a phage sequence and nothing here
says anything about phage biology.
"""

import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_benchmark import (  # noqa: E402
    concat,
    evaluate,
    read_fasta,
    revcomp,
    run,
    wgsim,
    write_fasta,
)

# A tRNA gene is about 76 bp, but integration disrupts a region around it and
# the reads that matter are the ones spanning the junction, so the repeated
# unit here is the length of a read pair's footprint rather than of a tRNA.
REPEAT = 300
# Real tRNA families in E. coli run to roughly 80 copies.
COPIES = 80


def build_host(host_seq, rng, copies=COPIES):
    """Put `copies` near-identical repeats into the host at spread positions.

    Returns the new sequence and the position of each repeat in it. The copies
    differ at a few bases, as members of a real family do, which is enough to
    stop them being literally identical without making them distinguishable to
    an aligner.
    """
    unit = "".join(rng.choice("ACGT") for _ in range(REPEAT))
    spacing = len(host_seq) // (copies + 1)
    pieces, positions, cursor = [], [], 0
    for i in range(1, copies + 1):
        at = i * spacing
        copy = list(unit)
        for _ in range(3):
            j = rng.randrange(REPEAT)
            copy[j] = rng.choice([b for b in "ACGT" if b != copy[j]])
        pieces += [host_seq[cursor:at], "".join(copy)]
        positions.append(sum(len(p) for p in pieces))
        cursor = at
    pieces.append(host_seq[cursor:])
    return "".join(pieces), positions


def build_sample(host_seq, repeat_positions, element_seq, args, rng):
    """Insert element fragments, either just after a repeat or in unique sequence."""
    if args.where == "adjacent":
        chosen = sorted(rng.sample(repeat_positions, args.n_insertions))
    else:
        # Unique sequence, kept well clear of every repeat.
        chosen = []
        while len(chosen) < args.n_insertions:
            pos = rng.randint(10_000, len(host_seq) - 10_000)
            if min(abs(pos - r) for r in repeat_positions) < 5_000:
                continue
            if all(abs(pos - p) > 20_000 for p in chosen):
                chosen.append(pos)
        chosen.sort()

    pieces, truth, cursor = [], [], 0
    for pos in chosen:
        length = rng.randint(args.insert_min, args.insert_max)
        start = rng.randint(0, len(element_seq) - length)
        fragment = element_seq[start : start + length]
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        pieces += [host_seq[cursor:pos], fragment]
        truth.append(("host", pos, length))
        cursor = pos
    pieces.append(host_seq[cursor:])
    return "".join(pieces), truth


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--host", required=True, help="Bacterial host genome FASTA")
    p.add_argument("--element", required=True, help="Genome the inserted element is cut from")
    p.add_argument("--contaminants", nargs="+", required=True)
    p.add_argument("--kraken-db", required=True)
    p.add_argument("--nodes", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--host-taxid", type=int, default=562, help="562 is Escherichia coli")
    p.add_argument("--where", choices=("adjacent", "unique"), default="adjacent")
    p.add_argument("--n-insertions", type=int, default=20)
    p.add_argument("--insert-min", type=int, default=2_000)
    p.add_argument("--insert-max", type=int, default=8_000)
    p.add_argument("--coverage", type=float, default=30.0)
    p.add_argument("--contaminant-pairs", type=int, default=200_000)
    p.add_argument("--chimera-pairs", type=int, default=10_000)
    p.add_argument("--repeat-copies", type=int, default=COPIES)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    (_, raw_host), = read_fasta(args.host).items()
    element_seq = max(read_fasta(args.element).values(), key=len)
    if len(element_seq) < args.insert_max:
        raise SystemExit(
            f"The element is {len(element_seq)} bp but insertions are up to {args.insert_max} bp."
        )

    host_seq, repeat_positions = build_host(raw_host.upper(), rng, args.repeat_copies)
    if args.n_insertions > len(repeat_positions):
        raise SystemExit(
            f"{args.n_insertions} insertions asked for but only {len(repeat_positions)} "
            f"repeat copies exist; raise --repeat-copies or lower --n-insertions."
        )
    reference = out / "host.fa"
    write_fasta(reference, {"host": host_seq})
    print(f"host {len(host_seq):,} bp with {len(repeat_positions)} copies of a {REPEAT} bp repeat")

    sample_seq, truth = build_sample(host_seq, repeat_positions, element_seq, args, rng)
    sample = out / "sample.fa"
    write_fasta(sample, {"sample": sample_seq})
    with open(out / "truth.bed", "w") as handle:
        for chrom, pos, length in truth:
            handle.write(f"{chrom}\t{pos}\t{pos + 1}\t{length}\n")
    print(f"{len(truth)} insertions, {args.where} to a repeat copy")

    # Reads: the sample, the contaminants, and the chance chimeras.
    pairs = int(args.coverage * len(sample_seq) / 300)
    wgsim(sample, pairs, out / "s1.fq", out / "s2.fq", args.seed)
    contaminants = {n: s for path in args.contaminants for n, s in read_fasta(path).items()}
    per = args.contaminant_pairs // max(1, len(args.contaminants))
    reads1, reads2 = [out / "s1.fq"], [out / "s2.fq"]
    for i, path in enumerate(args.contaminants):
        wgsim(path, per, out / f"c{i}_1.fq", out / f"c{i}_2.fq", args.seed + 100 + i)
        reads1.append(out / f"c{i}_1.fq")
        reads2.append(out / f"c{i}_2.fq")
    chim1, chim2 = out / "ch1.fq", out / "ch2.fq"
    _write_chimeras(host_seq, contaminants, args.chimera_pairs, rng, chim1, chim2)
    reads1.append(chim1)
    reads2.append(chim2)
    concat(reads1, out / "all_1.fq")
    concat(reads2, out / "all_2.fq")

    run(["bwa", "index", str(reference)])
    with open(out / "aligned.sam", "w") as handle:
        subprocess.run(
            ["bwa", "mem", "-t", str(args.threads), str(reference),
             str(out / "all_1.fq"), str(out / "all_2.fq")],
            stdout=handle, stderr=subprocess.DEVNULL, check=True,
        )
    run(["samtools", "sort", "-o", str(out / "aligned.bam"), str(out / "aligned.sam")])
    run(["samtools", "index", str(out / "aligned.bam")])

    run(["kraken2", "--db", args.kraken_db, "--threads", str(args.threads),
         "--paired", "--output", str(out / "kraken2.txt"),
         "--report", str(out / "kraken2.report"),
         str(out / "all_1.fq"), str(out / "all_2.fq")])

    run(["hgt-site", "call", "--bam", str(out / "aligned.bam"),
         "--classifications", str(out / "kraken2.txt"),
         "--host-taxid", str(args.host_taxid), "--nodes", args.nodes,
         "--threads", str(args.threads), "-o", str(out / "sites.tsv")])

    import csv
    with open(out / "sites.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    calls = [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows if r["filter"] == "PASS"]
    stats = evaluate(calls, [(c, p, 0) for c, p, _ in truth])

    from collections import Counter
    flags = Counter(r["filter"] for r in rows)
    print(f"\ncondition          {args.where}")
    print(f"candidates         {len(rows)}")
    print(f"sites reported     {stats['calls']}")
    print(f"of those real      {stats['true_calls']}")
    print(f"precision          {stats['precision']:.3f}")
    print(f"recall             {stats['recall']:.3f} ({stats['recovered']} of {stats['truth']})")
    print("\nfilter column:")
    for flag, n in flags.most_common(8):
        print(f"  {flag:<28} {n}")


def _write_chimeras(host_seq, contaminants, n, rng, out1, out2, molecule=500, length=150):
    """Host-plus-contaminant chimeric pairs, the artefact the test is about."""
    names = list(contaminants)
    with open(out1, "w") as h1, open(out2, "w") as h2:
        for i in range(n):
            half = molecule // 2
            hpos = rng.randint(0, len(host_seq) - half - 1)
            host_half = host_seq[hpos : hpos + half]
            source = contaminants[rng.choice(names)]
            cpos = rng.randint(0, len(source) - half - 1)
            other_half = source[cpos : cpos + half]
            molecule_seq = host_half + other_half
            read1 = molecule_seq[:length]
            read2 = revcomp(molecule_seq[-length:])
            for handle, seq, mate in ((h1, read1, 1), (h2, read2, 2)):
                handle.write(f"@chimera_{i}/{mate}\n{seq}\n+\n{'I' * len(seq)}\n")


if __name__ == "__main__":
    main()
