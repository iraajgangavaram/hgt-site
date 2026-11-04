# Installation

## From PyPI

```bash
conda install -c conda-forge -c bioconda hgt-site
```

or from PyPI:

```bash
pip install hgt-site
```

The Bioconda package comes with a container built automatically, so hgt-site can be run with nothing installed:

```bash
docker run --rm quay.io/biocontainers/hgt-site:0.3.0--pyhdfd78af_0 hgt-site --help
```

The exact tag for a release is listed on its [Bioconda page](https://bioconda.github.io/recipes/hgt-site/README.html), and the matching Singularity image is on the [Galaxy depot](https://depot.galaxyproject.org/singularity/).

## From Bioconda

```bash
conda install -c conda-forge -c bioconda hgt-site
```

Bioconda packages are built and tested on their own infrastructure before they
are published, so a successful install is also a check that the package works
on a clean machine.

## From source

```bash
git clone https://github.com/iraajgangavaram/hgt-site.git
cd hgt-site
pip install .
```

## Requirements

- Python 3.10 or newer.
- pysam, installed as a dependency. It is needed for BAM and CRAM input; plain SAM files can be read without it.

## Checking the install

```bash
hgt-site --version
```

## Check it works

```bash
hgt-site selftest
```

hgt-site builds a small dataset with known insertion sites, calls it, and checks that what comes back is what went in. It needs no data, no Kraken2 database and no network, and it leaves nothing behind. If it prints that hgt-site is working, the installation is sound and any later problem is in the inputs rather than the install.
