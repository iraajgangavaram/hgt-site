"""Mutate a genome to a target nucleotide divergence.

Used by the divergence benchmark: the sample's reads are simulated from the
mutated copy while the Kraken2 database stays built from the original, so the
classifier is looking for a strain it has never seen.

    python benchmark/mutate_genome.py ecoli.fa ecoli_5pc.fa 0.05 11

The four arguments are the source FASTA, where to write it, the fraction of
bases to substitute, and a seed. Substitutions are drawn uniformly across the
sequence, which is not what real divergence looks like: real divergence is
clustered, carries indels, and includes whole genes absent from any reference.
It is the harder case for a k-mer classifier in one respect, since there are
no long untouched stretches, and the easier case in another, since no sequence
is missing altogether.
"""

import random
import sys

BASES = "ACGT"


def read_fasta(path):
    name, seq, out = None, [], {}
    for line in open(path):
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name:
        out[name] = "".join(seq)
    return out


def mutate(seq, rate, rng):
    """Substitute `rate` of the bases, never for themselves."""
    out = list(seq.upper())
    for i in rng.sample(range(len(out)), int(len(out) * rate)):
        if out[i] in BASES:
            out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def main(argv):
    if len(argv) != 5:
        print(__doc__.strip().split("\n\n")[2], file=sys.stderr)
        return 2
    source, destination, rate, seed = argv[1], argv[2], float(argv[3]), int(argv[4])
    if not 0 <= rate <= 1:
        print(f"divergence must be between 0 and 1, not {rate}", file=sys.stderr)
        return 2
    rng = random.Random(seed)
    with open(destination, "w") as handle:
        for name, seq in read_fasta(source).items():
            handle.write(f">{name}\n")
            mutated = mutate(seq, rate, rng)
            for i in range(0, len(mutated), 80):
                handle.write(mutated[i : i + 80] + "\n")
    print(f"{destination}: {rate * 100:g}% divergence from {source}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
