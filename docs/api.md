# Python API

Everything the command line does is available as a library, so hgt-site can be
used from a notebook or inside another program without shelling out.

```python
from hgt_site import call_sites, passed, to_dataframe

sites = call_sites("sample.bam", ["classified.txt"], host_taxid=9606)
good = passed(sites)
print(f"{len(sites)} candidates, {len(good)} passed")

frame = to_dataframe(good)
frame[["chrom", "start", "n_pairs", "expected_pairs", "qvalue"]]
```

`to_dataframe` needs pandas, which hgt-site does not install for you. Nothing
else in the API does.

## call_sites

```python
call_sites(
    alignments,            # BAM, CRAM or SAM path
    kraken,                # list of per-read classification files
    host_taxid,            # 9606 for human
    nodes_dmp=None,        # NCBI nodes.dmp, to exclude the whole host lineage
    exclude=(),            # extra taxids to ignore
    reference=None,        # reference FASTA, needed for CRAM
    threads=1,
    split_reads=True,      # find exact junctions in soft-clipped reads
    annotation=None,       # GTF or GFF3, to name the gene at each site
    max_gene_distance=100_000,
    classifier="auto",     # kraken2, kaiju, centrifuge, diamond, blast, tsv
    summary=None,          # pass a dict to have read counts added to it
    **thresholds,          # any key of DEFAULTS
)
```

Returns a list of `Site` objects, every candidate rather than only the ones
that passed, so a rejected position can be inspected. Each has `chrom`,
`start`, `end`, a `flags` list and a `stats` dict holding every column of the
output table.

The thresholds are the same ones the command line exposes:

```python
from hgt_site import DEFAULTS

DEFAULTS            # {'min_mapq': 20, 'max_gap': 300, ...}
call_sites(..., min_support=5, max_qvalue=0.01)
```

A `summary` dict is filled in with the read counts behind the run, which is
what the MultiQC output reports:

```python
summary = {}
sites = call_sites("sample.bam", ["classified.txt"], 9606, summary=summary)
summary["anchored_pairs"], summary["foreign_anchored_pairs"]
```

## Writing the output formats

```python
from hgt_site import write_table, write_bed, write_vcf, write_multiqc

write_table(sites, "sites.tsv")
write_bed(sites, "sites.bed")
write_vcf(sites, "sites.vcf.gz", version="0.3.0", contigs={"chr21": 46709983})
write_multiqc(sites, "sample_mqc.json", sample="sample01")
```

## Reading classifications on their own

```python
from hgt_site import load_classifications

calls = load_classifications(["classified.txt"])     # {read name: taxid}
calls = load_classifications(["hits.tsv"], "blast")  # when detection gets it wrong
```

## Making a test dataset

```python
from hgt_site import simulate
from hgt_site.simulate import SimConfig

truth = simulate("demo", SimConfig(seed=1))
```

That writes `demo.sam`, `demo.kraken2.txt`, `demo.truth.bed` and `demo.host.fa`,
and returns the insertions it put in, so a change can be checked against known
answers. This is what `hgt-site selftest` runs.
