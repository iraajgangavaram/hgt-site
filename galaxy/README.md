# Galaxy wrapper

A Galaxy tool wrapper for `hgt-site call`, so the tool can be run by people who
do not use a command line. Galaxy installs the underlying program from Bioconda,
which is why this could not be submitted until the recipe had merged.

```
tools/hgt_site/
    hgt_site.xml          the tool definition
    .shed.yml             Tool Shed metadata
    test-data/            a small simulated dataset with known insertion sites
```

## Testing it

With [planemo](https://planemo.readthedocs.io/):

```bash
planemo lint galaxy/tools/hgt_site
planemo test --biocontainers galaxy/tools/hgt_site
planemo serve --biocontainers galaxy/tools/hgt_site   # to try it in a browser
```

`planemo test` installs hgt-site from Bioconda, runs the tool on the dataset in
`test-data/`, and checks the outputs. The test data is a simulated sample with
two detectable insertions from *Escherichia coli*, both of which the tool should
report, which is why the BED assertion expects exactly two lines.

## Submitting it

Galaxy tools are distributed through the Tool Shed, and the usual route is a
pull request to [galaxyproject/tools-iuc](https://github.com/galaxyproject/tools-iuc),
which maintains tools collectively rather than leaving each one to its author.
Copy `tools/hgt_site/` into `tools/` in a fork of that repository, check
`planemo lint` is clean, and open a pull request. Their CI runs the tests on
every pull request, so a wrapper that passes locally usually passes there.

If the IUC decline it, the same files can go to a personal Tool Shed account
instead, with `owner` in `.shed.yml` changed from `iuc` to that account.

## Keeping it up to date

The wrapper pins a version in one place:

```xml
<token name="@TOOL_VERSION@">0.3.0</token>
```

Galaxy requires that version to exist on Bioconda, so the order for a release is
always: tag the release, update the Bioconda recipe, wait for it to merge and
build, then bump this token here.

The `+galaxyN` suffix on the tool version is separate. It goes back to
`+galaxy0` whenever `@TOOL_VERSION@` changes, and increments to `+galaxy1`,
`+galaxy2` and so on when the wrapper itself changes while the underlying
version stays the same.
