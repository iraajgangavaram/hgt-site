"""Do exon junctions turn into false insertion sites?

The companion script, rnaseq_test.py, samples reads from contiguous genomic
sequence, so nothing in it crosses an exon junction. That leaves the question
that matters most for RNA-seq untested: hgt-site reads soft clips to place a
junction to the base pair, and a read crossing an exon junction aligned with
a DNA aligner comes back soft-clipped at exactly that boundary.

This builds genes with introns, samples reads from the spliced transcript,
and aligns them to the genome with bwa, which knows nothing about splicing.
Every transcript read that crosses a junction is therefore soft-clipped, with
the clipped part being real host sequence from the next exon. That is the
pessimistic case: anyone using STAR or HISAT2 gets a cleaner alignment than
this, so a result here is a floor rather than an estimate.

Two things are being looked for:

  - False sites at exon boundaries. A soft clip only matters if the position
    already has enriched foreign support, so in principle an exon boundary
    with no foreign reads should never become a candidate. This checks that
    in principle holds in practice.
  - Junctions placed at an exon boundary rather than at the insertion, for
    sites that are real.
"""

import argparse
import csv
import json
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_benchmark import evaluate, read_fasta, revcomp, run, write_fasta  # noqa: E402

READ_LEN = 150
FRAGMENT = 400
ERROR_RATE = 0.002
HOST_TAXID = 9606


def log_uniform(rng, low, high):
    import math

    return math.exp(rng.uniform(math.log(low), math.log(high)))


def build_genes(host_len, n_genes, max_fold, rng):
    """Genes of three to five exons, with introns between them."""
    spacing = host_len // n_genes
    genes = []
    for i in range(n_genes):
        base = i * spacing
        exons, cursor = [], base + rng.randint(0, 200)
        for _ in range(rng.randint(3, 5)):
            length = rng.randint(400, 1_500)
            exons.append((cursor, length))
            cursor += length + rng.randint(500, 2_000)
        if cursor - base > spacing:
            continue
        genes.append({"exons": exons, "expression": log_uniform(rng, 1, max_fold)})
    return genes


def transcript(host_seq, gene):
    """The spliced sequence, and the offset of each exon boundary within it."""
    pieces, boundaries, total = [], [], 0
    for start, length in gene["exons"]:
        pieces.append(host_seq[start : start + length])
        total += length
        boundaries.append(total)
    return "".join(pieces), boundaries[:-1]


def mutate(seq, rng):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < ERROR_RATE:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def emit_pair(name, seq, rng, f1, f2):
    if len(seq) < READ_LEN:
        return
    fragment = min(FRAGMENT, len(seq))
    offset = rng.randint(0, len(seq) - fragment) if len(seq) > fragment else 0
    piece = seq[offset : offset + fragment]
    qual = "I" * READ_LEN
    f1.write(f"@{name}/1\n{mutate(piece[:READ_LEN], rng)}\n+\n{qual}\n")
    f2.write(f"@{name}/2\n{mutate(revcomp(piece)[:READ_LEN], rng)}\n+\n{qual}\n")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--host", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--contaminants", nargs="+", required=True)
    p.add_argument("--kraken-db", required=True)
    p.add_argument("--nodes", required=True)
    p.add_argument("--genes", type=int, default=400)
    p.add_argument("--max-fold", type=float, default=10_000)
    p.add_argument("--pairs", type=int, default=300_000)
    p.add_argument("--insertions", type=int, default=20)
    p.add_argument("--insert-support", type=int, default=12)
    p.add_argument("--chimera-pairs", type=int, default=10_000)
    p.add_argument("--donor-fragment", type=int, default=5_000)
    p.add_argument("--local-background", action="store_true")
    p.add_argument("--exclude-exon-boundaries", type=int, default=0,
                   help="write a BED of every exon boundary with this much slop either side and pass "
                        "it to --exclude-regions, which is what a user could build from a GTF")
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
    expressed = sorted(genes, key=lambda g: -g["expression"])[: max(len(genes) // 2, args.insertions)]
    chosen = rng.sample(expressed, args.insertions)
    truth, exon_boundaries = [], []
    for gene in chosen:
        # The insertion goes inside the middle exon, away from its edges.
        exon = gene["exons"][len(gene["exons"]) // 2]
        # At least 200 bp of exon on the host side, so the anchor can map.
        offset = rng.randint(200, exon[1] - 200)
        gene["insertion"] = (exon[0], offset)
        truth.append((chrom, exon[0] + offset, gene["expression"]))

    for gene in genes:
        for start, length in gene["exons"][:-1]:
            exon_boundaries.append(start + length)

    weights = [g["expression"] * sum(length for _, length in g["exons"]) for g in genes]
    total = sum(weights)

    with open(out / "reads_1.fq", "w") as f1, open(out / "reads_2.fq", "w") as f2:
        for i, gene in enumerate(genes):
            seq, _ = transcript(host_seq, gene)
            if "insertion" in gene:
                exon_start, offset = gene["insertion"]
                # Where the insertion sits within the spliced transcript.
                before = 0
                for start, length in gene["exons"]:
                    if start == exon_start:
                        break
                    before += length
                cut = before + offset
                seq = seq[:cut] + donor_seq[: args.donor_fragment] + seq[cut:]
            n = int(args.pairs * weights[i] / total)
            for j in range(n):
                emit_pair(f"g{i}_{j}", seq, rng, f1, f2)
        for k, gene in enumerate(chosen):
            exon_start, offset = gene["insertion"]
            left = host_seq[exon_start : exon_start + offset]
            for j in range(args.insert_support):
                host_part = left[-rng.randint(150, min(200, max(len(left), 151))):]
                emit_pair(f"j{k}_{j}", host_part + donor_seq[: 500 - len(host_part)], rng, f1, f2)
        for i in range(args.chimera_pairs):
            gene = rng.choices(genes, weights=weights, k=1)[0]
            seq, _ = transcript(host_seq, gene)
            if len(seq) < 320:
                continue
            hstart = rng.randint(0, len(seq) - 300)
            host_part = seq[hstart : hstart + rng.randint(250, 300)]
            source = rng.choice(contam)
            cstart = rng.randint(0, len(source) - 500)
            molecule = host_part + source[cstart : cstart + 500 - len(host_part)]
            emit_pair(f"chim{i}", revcomp(molecule) if rng.random() < 0.5 else molecule, rng, f1, f2)

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
    if args.exclude_exon_boundaries:
        slop = args.exclude_exon_boundaries
        bed = out / "exon_boundaries.bed"
        with open(bed, "w") as fh:
            for b in sorted(exon_boundaries):
                fh.write(f"{chrom}\t{max(b - slop, 0)}\t{b + slop}\n")
        cmd += ["--exclude-regions", str(bed)]
    run(cmd)

    with open(out / "sites.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    passed = [r for r in rows if r["filter"] == "PASS"]
    calls = [(r["chrom"], int(r["start"]), int(r["end"])) for r in passed]
    result = evaluate(calls, truth)

    # How many of the sites that are not real sit on an exon boundary, and
    # how many real sites had their junction pulled to one.
    bset = set(exon_boundaries)
    truth_pos = {p for _, p, _ in truth}
    false_on_boundary = 0
    junction_on_boundary = 0
    for r in passed:
        start, end = int(r["start"]), int(r["end"])
        is_true = any(start - 500 <= p <= end + 500 for p in truth_pos)
        near_boundary = any(start - 200 <= b <= end + 200 for b in bset)
        if not is_true and near_boundary:
            false_on_boundary += 1
        for key in ("left_junction", "right_junction"):
            junction = r.get(key) or ""
            if junction.isdigit() and int(junction) in bset:
                junction_on_boundary += 1
                break

    summary = {
        "seed": args.seed,
        "spliced": True,
        "genes": len(genes),
        "exon_boundaries": len(exon_boundaries),
        "local_background": args.local_background,
        "exclude_exon_boundaries": args.exclude_exon_boundaries,
        "candidates": len(rows),
        "hgt_site": result,
        "false_sites_on_an_exon_boundary": false_on_boundary,
        "junctions_placed_on_an_exon_boundary": junction_on_boundary,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in
                      ("hgt_site", "false_sites_on_an_exon_boundary",
                       "junctions_placed_on_an_exon_boundary", "exon_boundaries")}))


if __name__ == "__main__":
    main()
