"""Does hgt-site work on RNA-seq once the background follows coverage?

RNA-seq breaks the even-coverage assumption far harder than a copy-number
gain does. Coverage follows expression across several orders of magnitude,
there is no coverage at all between genes, and the chimeras made during
library preparation come from whatever transcripts are abundant. If the
even-coverage background is the only thing standing in the way, then
--local-background should rescue it; if something else is wrong, this is
where it shows.

What this simulates:
  - A host region split into genes, each with an expression level drawn
    log-uniformly across four orders of magnitude.
  - Reads sampled from genes in proportion to expression times length, with
    nothing between genes, so total depth is a fixed budget as it is in a
    real library.
  - Insertions placed inside genes, since an insertion in a gene nobody
    transcribes produces no reads and cannot be found by any method.
  - Chance chimeras whose host half is drawn from the same expression
    profile, because an abundant transcript makes abundant chimeras.

What it does not simulate: splicing. Reads are sampled from contiguous
genomic sequence, so no read crosses an exon junction. A real spliced
aligner, and a real intron, are not tested here, and that limitation is the
main reason this is a feasibility check rather than a benchmark.
"""

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_benchmark import evaluate, read_fasta, revcomp, run, write_fasta  # noqa: E402

GENE_LEN = 2_000
READ_LEN = 150
FRAGMENT = 400
ERROR_RATE = 0.002
DONOR_TAXID = 562  # Escherichia coli, as in the other benchmarks
HOST_TAXID = 9606


def log_uniform(rng, low, high):
    import math

    return math.exp(rng.uniform(math.log(low), math.log(high)))


def build_genes(host_len, n_genes, max_fold, rng):
    """Evenly spaced genes with log-uniform expression."""
    spacing = host_len // n_genes
    genes = []
    for i in range(n_genes):
        start = i * spacing + rng.randint(0, max(spacing - GENE_LEN - 1, 1))
        genes.append({"start": start, "length": GENE_LEN, "expression": log_uniform(rng, 1, max_fold)})
    return genes


def mutate(seq, rng):
    if ERROR_RATE <= 0:
        return seq
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < ERROR_RATE:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def emit_pair(name, seq, start, rng, fastq1, fastq2):
    """One read pair from seq at a random fragment, with errors."""
    if len(seq) < READ_LEN:
        return
    fragment = min(FRAGMENT, len(seq))
    offset = rng.randint(0, len(seq) - fragment) if len(seq) > fragment else 0
    piece = seq[offset : offset + fragment]
    r1 = mutate(piece[:READ_LEN], rng)
    r2 = mutate(revcomp(piece)[:READ_LEN], rng)
    qual = "I" * READ_LEN
    fastq1.write(f"@{name}/1\n{r1}\n+\n{qual}\n")
    fastq2.write(f"@{name}/2\n{r2}\n+\n{qual}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--host", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--contaminants", nargs="+", required=True)
    p.add_argument("--kraken-db", required=True)
    p.add_argument("--nodes", required=True)
    p.add_argument("--genes", type=int, default=400)
    p.add_argument("--max-fold", type=float, default=10_000)
    p.add_argument("--pairs", type=int, default=400_000, help="host read pairs in total")
    p.add_argument("--insertions", type=int, default=20)
    p.add_argument("--insert-support", type=int, default=12, help="pairs across each junction")
    p.add_argument("--chimera-pairs", type=int, default=10_000)
    p.add_argument("--donor-fragment", type=int, default=5_000)
    p.add_argument("--local-background", action="store_true")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    host = read_fasta(args.host)
    chrom, host_seq = next(iter(host.items()))
    donor_seq = next(iter(read_fasta(args.donor).values()))
    contam = [s for g in args.contaminants for s in read_fasta(g).values() if len(s) > 10_000]

    genes = build_genes(len(host_seq), args.genes, args.max_fold, rng)
    # Insertions go in genes that are actually transcribed. A gene in the
    # bottom half of the expression range would yield almost no reads, so
    # there would be nothing for any method to find.
    expressed = sorted(genes, key=lambda g: -g["expression"])[: max(args.genes // 2, args.insertions)]
    chosen = rng.sample(expressed, args.insertions)
    truth = []
    for gene in chosen:
        offset = rng.randint(400, gene["length"] - 400)
        gene["insertion"] = offset
        truth.append((chrom, gene["start"] + offset, gene["expression"]))

    weights = [g["expression"] * g["length"] for g in genes]
    total_weight = sum(weights)

    with open(out / "reads_1.fq", "w") as f1, open(out / "reads_2.fq", "w") as f2:
        # Host reads, budget split by expression.
        for i, gene in enumerate(genes):
            n = int(args.pairs * weights[i] / total_weight)
            seq = host_seq[gene["start"] : gene["start"] + gene["length"]]
            if "insertion" in gene:
                cut = gene["insertion"]
                seq = seq[:cut] + donor_seq[: args.donor_fragment] + seq[cut:]
            for j in range(n):
                emit_pair(f"g{i}_{j}", seq, gene["start"], rng, f1, f2)
        # Extra pairs right across each junction, so a real insertion has
        # support whatever its expression level happens to be.
        for k, gene in enumerate(chosen):
            cut = gene["insertion"]
            left = host_seq[gene["start"] : gene["start"] + cut]
            for j in range(args.insert_support):
                host_part = left[-rng.randint(120, 200):]
                donor_part = donor_seq[: 500 - len(host_part)]
                emit_pair(f"j{k}_{j}", host_part + donor_part, 0, rng, f1, f2)
        # Chance chimeras, host half drawn from the same expression profile.
        for i in range(args.chimera_pairs):
            gene = rng.choices(genes, weights=weights, k=1)[0]
            hstart = gene["start"] + rng.randint(0, gene["length"] - 300)
            host_part = host_seq[hstart : hstart + rng.randint(250, 300)]
            source = rng.choice(contam)
            cstart = rng.randint(0, len(source) - 500)
            seq = host_part + source[cstart : cstart + 500 - len(host_part)]
            if rng.random() < 0.5:
                seq = revcomp(seq)
            emit_pair(f"chim{i}", seq, 0, rng, f1, f2)

    write_fasta(out / "host.fa", {chrom: host_seq})
    run(["bwa", "index", out / "host.fa"])
    with open(out / "aligned.sam", "w") as sam:
        subprocess.run(
            ["bwa", "mem", "-t", str(args.threads), str(out / "host.fa"),
             str(out / "reads_1.fq"), str(out / "reads_2.fq")],
            stdout=sam, stderr=subprocess.DEVNULL, check=True,
        )
    run(["samtools", "sort", "-@", str(args.threads), "-o", out / "full.bam", out / "aligned.sam"])
    run(["samtools", "index", out / "full.bam"])
    # Unmapped mates, for the classifier.
    with open(out / "unmapped.fq", "w") as fq:
        subprocess.run(["samtools", "fastq", "-f", "4", str(out / "full.bam")],
                       stdout=fq, stderr=subprocess.DEVNULL, check=True)
    run(["kraken2", "--db", args.kraken_db, "--threads", str(args.threads),
         "--output", out / "kraken2.txt", "--report", out / "kraken2.report",
         str(out / "unmapped.fq")])

    cmd = [sys.executable, "-m", "hgt_site", "call", "--bam", str(out / "full.bam"),
           "--classifications", str(out / "kraken2.txt"), "--host-taxid", str(HOST_TAXID),
           "--nodes", args.nodes, "-o", str(out / "sites.tsv")]
    if args.local_background:
        cmd.append("--local-background")
    run(cmd)

    import csv

    with open(out / "sites.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    calls = [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows if r["filter"] == "PASS"]
    result = evaluate(calls, truth)
    summary = {
        "seed": args.seed,
        "genes": args.genes,
        "max_fold": args.max_fold,
        "pairs": args.pairs,
        "chimera_pairs": args.chimera_pairs,
        "local_background": args.local_background,
        "hgt_site": result,
        "candidates": len(rows),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["hgt_site"]))


if __name__ == "__main__":
    main()
