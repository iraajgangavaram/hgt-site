"""Benchmark hgt-site on reads simulated from real genomes.

Steps:
1. Insert fragments of a donor genome into a host sequence at known positions.
2. Simulate paired reads from that sample with wgsim.
3. Add free contaminant read pairs and host-contaminant chimeric pairs.
4. Align to the host with BWA-MEM, classify unmapped reads with Kraken2.
5. Run hgt-site and compare its calls with the known insertion sites.

Needs bwa, samtools, wgsim and kraken2 on PATH, a BWA index of the host and
a Kraken2 database that contains the host and the bacterial genomes.
"""

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    records, name, chunks = {}, None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    records[name] = "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name:
        records[name] = "".join(chunks).upper()
    return records


def write_fasta(path, records):
    with open(path, "w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i : i + 60] + "\n")


def run(cmd, **kwargs):
    if not isinstance(cmd, str):
        cmd = [str(c) for c in cmd]
    print("+", cmd if isinstance(cmd, str) else " ".join(map(str, cmd)), file=sys.stderr)
    subprocess.run(cmd, check=True, **kwargs)


def build_sample(host, donor, args, rng):
    """Return the sample sequence and the true insertion sites in host coordinates."""
    (host_name, host_seq), = host.items()
    donor_seq = max(donor.values(), key=len)
    positions = []
    while len(positions) < args.n_insertions:
        pos = rng.randint(10_000, len(host_seq) - 10_000)
        if all(abs(pos - p) > 20_000 for p in positions):
            positions.append(pos)
    positions.sort()
    pieces, truth, cursor = [], [], 0
    for pos in positions:
        length = rng.randint(args.insert_min, args.insert_max)
        start = rng.randint(0, len(donor_seq) - length)
        fragment = donor_seq[start : start + length]
        if rng.random() < 0.5:
            fragment = revcomp(fragment)
        pieces += [host_seq[cursor:pos], fragment]
        truth.append((host_name, pos, length))
        cursor = pos
    pieces.append(host_seq[cursor:])
    return host_name, "".join(pieces), truth


def build_sample_intvi(args, out):
    """Insert donor fragments with the published simulator of Scott et al. (2022).

    Uses scripts/insert_virus.py from aehrc/vector-integration-simulation-pipeline,
    which adds rearranged and deleted fragments, host deletions and junction gaps
    or overlaps. Whole-genome insertions are switched off, since the donor here is
    a bacterium, and no free episomes are added.
    """
    import csv as _csv

    run(
        [sys.executable, args.intvi_script, "--host", args.host, "--virus", args.donor,
         "--ints", out / "sample.fa", "--int_info", out / "intvi_ints.tsv",
         "--epi_num", "0", "--epi_info", out / "intvi_epi.tsv",
         "--int_num", args.n_insertions, "--p_whole", "0",
         "--min-len", args.insert_min, "--max-len", args.insert_max,
         "--min-sep", "20000", "--seed", args.seed],
        stdout=subprocess.DEVNULL,
    )
    sample = read_fasta(out / "sample.fa")
    (host_name,) = read_fasta(args.host).keys()
    truth = []
    with open(out / "intvi_ints.tsv") as handle:
        for row in _csv.DictReader(handle, delimiter="\t"):
            truth.append((row["chr"], int(row["hPos_input_fasta"]), 0))
    return host_name, sample[host_name], truth


def build_chimeras(host_seq, contaminants, n, rng, molecule=500):
    """Host-contaminant molecules of the kind made by chance during library preparation.

    Molecules are 500 bp because wgsim skips sequences shorter than the
    insert size plus three standard deviations (460 bp here).
    """
    contam = [s for genome in contaminants for s in genome.values() if len(s) > 10_000]
    molecules = {}
    for i in range(n):
        host_len = rng.randint(250, 300)
        hstart = rng.randint(0, len(host_seq) - host_len)
        source = rng.choice(contam)
        cstart = rng.randint(0, len(source) - molecule)
        seq = host_seq[hstart : hstart + host_len] + source[cstart : cstart + molecule - host_len]
        if rng.random() < 0.5:
            seq = revcomp(seq)
        molecules[f"chim{i}"] = seq
    return molecules


def wgsim(fasta, n_pairs, out1, out2, seed, insert=400, sd=40, length=150):
    run(
        ["wgsim", "-N", n_pairs, "-1", length, "-2", length, "-d", insert, "-s", sd,
         "-e", "0.002", "-r", "0.0005", "-R", "0", "-X", "0", "-S", seed, fasta, out1, out2],
        stdout=subprocess.DEVNULL,
    )


def concat(paths, out):
    with open(out, "w") as handle:
        for path in paths:
            handle.write(Path(path).read_text())


def evaluate(calls, truth, slop=500):
    matched, true_calls = set(), 0
    for chrom, start, end in calls:
        hits = [i for i, (c, p, _) in enumerate(truth) if c == chrom and start - slop <= p <= end + slop]
        if hits:
            true_calls += 1
            matched.update(hits)
    return {
        "calls": len(calls),
        "true_calls": true_calls,
        "truth": len(truth),
        "recovered": len(matched),
        "precision": true_calls / len(calls) if calls else 0.0,
        "recall": len(matched) / len(truth) if truth else 0.0,
    }


def lgtseek_calls(bam, min_pairs):
    """Group LGTSeek's host-side LGT reads into sites with hgt-site's clustering settings.

    LGTSeek reports read pairs, not sites, so the same clustering is applied to
    both tools and a site needs at least min_pairs pairs, as in the support-only
    baseline.
    """
    from hgt_site.alignments import REVERSE, Anchor
    from hgt_site.call import DEFAULTS
    from hgt_site.sites import cluster_anchors

    text = subprocess.run(["samtools", "view", str(bam)], capture_output=True, text=True, check=True).stdout
    anchors = []
    for line in text.splitlines():
        cols = line.split("\t")
        flag, start = int(cols[1]), int(cols[3]) - 1
        if flag & 0x4:
            continue
        length = sum(int(n) for n, op in __import__("re").findall(r"(\d+)([MDN=X])", cols[5]))
        anchors.append(Anchor(cols[0], cols[2], start, start + length, bool(flag & REVERSE), int(cols[4])))
    if not anchors:
        return []
    sites = cluster_anchors(anchors, DEFAULTS["max_gap"], DEFAULTS["max_span"])
    return [(s.chrom, s.start, s.end) for s in sites if len(s.anchors) >= min_pairs]


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--host", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--contaminants", nargs="+", required=True)
    p.add_argument("--kraken-db", required=True)
    p.add_argument("--host-taxid", type=int, default=9606)
    p.add_argument("--nodes", help="nodes.dmp passed on to hgt-site")
    p.add_argument("--outdir", required=True)
    p.add_argument("--n-insertions", type=int, default=20)
    p.add_argument("--insert-min", type=int, default=1000)
    p.add_argument("--insert-max", type=int, default=5000)
    p.add_argument("--coverage", type=float, default=10.0)
    p.add_argument("--contaminant-pairs", type=int, default=50_000)
    p.add_argument("--chimera-pairs", type=int, default=2_000)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--intvi-script", help="insert_virus.py from the Scott et al. simulator; used instead of the built-in insertion step")
    p.add_argument("--lgtseek-lib", help="LGTSeek perl5 lib directory; also runs LGTSeek's read classification")
    p.add_argument("--bacteria-index", help="BWA index of the bacterial reference genomes, for LGTSeek")
    p.add_argument("--min-lgtseek-pairs", type=int, default=3)
    p.add_argument("--keep-sorted-bam", action="store_true",
                   help="Keep a sorted, indexed copy of the host alignment (host.sorted.bam), for junctions.py")
    args = p.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    host = read_fasta(args.host)
    donor = read_fasta(args.donor)
    contaminants = [read_fasta(c) for c in args.contaminants]
    if args.intvi_script:
        host_name, sample_seq, truth = build_sample_intvi(args, out)
    else:
        host_name, sample_seq, truth = build_sample(host, donor, args, rng)
    write_fasta(out / "sample.fa", {host_name: sample_seq})
    with open(out / "truth.bed", "w") as handle:
        for chrom, pos, length in truth:
            handle.write(f"{chrom}\t{pos}\t{pos + 1}\tinsertion\t{length}\n")

    n_sample = int(args.coverage * len(sample_seq) / 300)
    wgsim(out / "sample.fa", n_sample, out / "s_1.fq", out / "s_2.fq", args.seed)

    contam_fa = out / "contaminants.fa"
    merged = {}
    for genome in contaminants:
        merged.update(genome)
    write_fasta(contam_fa, merged)
    wgsim(contam_fa, args.contaminant_pairs, out / "c_1.fq", out / "c_2.fq", args.seed + 1)

    (host_seq,) = host.values()
    write_fasta(out / "chimeras.fa", build_chimeras(host_seq, contaminants, args.chimera_pairs, rng))
    wgsim(out / "chimeras.fa", args.chimera_pairs, out / "x_1.fq", out / "x_2.fq", args.seed + 2, insert=400, sd=20)

    concat([out / "s_1.fq", out / "c_1.fq", out / "x_1.fq"], out / "reads_1.fq")
    concat([out / "s_2.fq", out / "c_2.fq", out / "x_2.fq"], out / "reads_2.fq")

    # Primary alignments only, in read order. LGTSeek needs this order, and
    # hgt-site ignores secondary and supplementary records anyway.
    run(
        f"bwa mem -t {args.threads} -v 1 {args.host} {out}/reads_1.fq {out}/reads_2.fq"
        f" | samtools view -b -F 0x900 -o {out}/host.bam -",
        shell=True,
    )
    # hgt-site only needs pairs with at least one unmapped mate, plus the header.
    run(f"samtools view -h -e 'flag.unmap || flag.munmap' -o {out}/aligned.sam {out}/host.bam", shell=True)
    run(f"samtools fastq -f 4 {out}/aligned.sam > {out}/unmapped.fq 2> /dev/null", shell=True)
    run(
        ["kraken2", "--db", args.kraken_db, "--threads", str(args.threads),
         "--output", str(out / "kraken2.txt"), "--report", str(out / "kraken2.report"), str(out / "unmapped.fq")],
        stderr=subprocess.DEVNULL,
    )

    cmd = [sys.executable, "-m", "hgt_site", "call", "--bam", str(out / "aligned.sam"),
           "--kraken", str(out / "kraken2.txt"), "--host-taxid", str(args.host_taxid),
           "-o", str(out / "sites.tsv")]
    if args.nodes:
        cmd += ["--nodes", args.nodes]
    run(cmd)

    import csv

    with open(out / "sites.tsv") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    passed = [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows if r["filter"] == "PASS"]
    support_only = [(r["chrom"], int(r["start"]), int(r["end"])) for r in rows if int(r["n_pairs"]) >= 3]
    if args.lgtseek_lib:
        run(
            f"bwa mem -t {args.threads} -v 1 {args.bacteria_index} {out}/reads_1.fq {out}/reads_2.fq"
            f" | samtools view -b -F 0x900 -o {out}/donor.bam -",
            shell=True,
        )
        script = Path(__file__).with_name("run_lgtseek.pl")
        samtools = subprocess.run(["which", "samtools"], capture_output=True, text=True).stdout.strip()
        run(["perl", script, args.lgtseek_lib, out / "donor.bam", out / "host.bam", out, samtools],
            stdout=open(out / "lgtseek_counts.txt", "w"), stderr=subprocess.DEVNULL)
        lgtseek_sites = lgtseek_calls(out / "lgtseek.lgt_recipient.bam", args.min_lgtseek_pairs)

    summary = {
        "seed": args.seed,
        "coverage": args.coverage,
        "contaminant_pairs": args.contaminant_pairs,
        "chimera_pairs": args.chimera_pairs,
        "hgt_site": evaluate(passed, truth),
        "support_only": evaluate(support_only, truth),
    }
    if args.lgtseek_lib:
        summary["lgtseek"] = evaluate(lgtseek_sites, truth)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))
    if args.keep_sorted_bam:
        run(["samtools", "sort", "-@", str(args.threads), "-o", out / "host.sorted.bam", out / "host.bam"])
        run(["samtools", "index", out / "host.sorted.bam"])
    for f in ["s_1.fq", "s_2.fq", "c_1.fq", "c_2.fq", "x_1.fq", "x_2.fq", "reads_1.fq", "reads_2.fq",
              "host.bam", "donor.bam"]:
        (out / f).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
