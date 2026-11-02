# Contributing

Bug reports, criticism of the method and runs on real data are all welcome.
There is no need to ask before opening an issue.

## The most useful thing you can send

A run on real sequencing data. hgt-site has been tested on simulated reads from
real genomes and not on real data with confirmed integrations, so a report from
a real sample is worth more than any code change. There is an issue template
for it, and a run that found nothing is as useful as one that worked.

## Reporting a problem

Use the issue templates. The INFO lines hgt-site prints are the most useful part
of a report, because they say how many anchored pairs it found and how many had
a Kraken2 call, which separates a problem with the input from a problem with the
calling.

If you disagree with the sites it reported, that is a bug report, not a
question.

## Changes to the code

```bash
git clone https://github.com/iraajgangavaram/hgt-site
cd hgt-site
pip install -e ".[test]"
pytest
ruff check src tests
```

The tests build small simulated datasets with known insertion sites and check
that the right sites come back, so a change that breaks the method fails the
tests rather than passing quietly. Please add a test with a change that alters
what gets called.

A few things to know before changing the calling itself:

- Results must not depend on the order of records in the file. A sorted BAM and
  a SAM of the same reads have to give the same sites, and there is a test for
  it. Anything that sorts or groups anchors needs a total ordering, not a
  partial one.
- Thresholds live in one place, `DEFAULTS` in `call.py`, and every one of them
  is exposed on the command line with its own help text.
- A new filter should add a flag rather than silently dropping a site. Every
  candidate stays in the table with the reason it failed, because a user
  checking a region needs to see that a site was considered and rejected.

## Changes to the benchmarks

The benchmark scripts are in `benchmark/`, and they rebuild their datasets from
seeds, so a result can be reproduced. If a change moves the numbers in the
README, the README has to move with it, including the confidence interval. A
claim in the README should say what was measured, on how many datasets, under
what conditions.

## Style

Plain British English in the documentation. Comments explain why something is
done, not what the line does. Lines up to 120 characters.

## Licence

Contributions are under the MIT licence, the same as the rest of the project.
