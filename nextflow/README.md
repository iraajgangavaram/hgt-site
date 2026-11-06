# Nextflow module

An nf-core module for `hgt-site call`, so the tool can be dropped into a
Nextflow pipeline without anyone having to find it first. It installs from
Bioconda, which is why it could not be written until the recipe merged.

```
modules/nf-core/hgtsite/call/
    main.nf            the process
    meta.yml           the module metadata
    environment.yml    the conda environment
    tests/main.nf.test an nf-test that chains Kraken2 to make its input
```

This is a standalone module and has nothing to do with any particular pipeline.
A pipeline that already classifies its unmapped reads can call it directly.

## Not yet tested

Neither Nextflow nor nf-test is available in the environment this was written
in, so **the module has not been run**. Before submitting it, run:

```bash
nf-core modules lint hgtsite/call
nf-test test modules/nf-core/hgtsite/call/tests/main.nf.test --profile docker --update-snapshot
```

from a checkout of [nf-core/modules](https://github.com/nf-core/modules) with
these files copied in. The snapshot file does not exist yet and is created by
that first run.

## About the test

nf-core/test-datasets has a Kraken2 database but no per-read classification
file, and a module test may not fetch data from anywhere else. So the test runs
`UNTAR` on the sarscov2 Kraken2 database, runs `KRAKEN2_KRAKEN2` on the reads
with `save_reads_assignment` set, and feeds that output to this module alongside
the matching alignment. That keeps everything inside the data nf-core already
ships.

The assertion is deliberately weak: it checks the table exists with its header,
rather than checking particular sites. The sarscov2 test data is not a host
sample with insertions in it, so no site should pass, and asserting on the
contents would be asserting on noise. The versions output is snapshotted.

## Submitting it

Fork [nf-core/modules](https://github.com/nf-core/modules), copy
`modules/nf-core/hgtsite/call/` into the same path there, make sure
`nf-core modules lint` is clean and the test passes, then open a pull request.
Their CI runs the test on every pull request.

## Keeping it up to date

Three places carry the version, and all three have to move together:

- `environment.yml`: `bioconda::hgt-site=0.3.0`
- `main.nf`: both container URLs

The version must exist on Bioconda first. The current pin is 0.3.0; it moves to
0.4.0 once that recipe merges and the container is built.

## Conventions this follows

nf-core changed how modules report versions. They no longer write a
`versions.yml` file; they emit a tuple on the `versions` topic with an `eval`
expression that is run by Nextflow:

```groovy
tuple val("${task.process}"), val('hgtsite'), eval('hgt-site --version | sed "s/^hgt-site //"'), topic: versions, emit: versions
```

If you are comparing this with an older module, that is why it looks different.
