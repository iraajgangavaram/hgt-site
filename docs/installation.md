# Installation

## From PyPI

```bash
pip install hgt-site
```

## From Bioconda

A recipe is under review. Once it is merged:

```bash
conda install -c conda-forge -c bioconda hgt-site
```

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
