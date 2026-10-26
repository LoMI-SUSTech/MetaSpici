# Validation — MetaSpici V1.1

## Checks performed for this release package

- CLI/package source updated to version `1.1`.
- Reference-mode workflow source retained from the user-validated reference baseline; the scientific defaults remain `--min-depth 0.15`, strict `FE > 50`, and `--min-contacts 4`.
- Short reads + Hi-C was reported by the maintainer as successfully running on their server.
- The MetaTOR environment was checked on the maintainer's server: `metator -h`, `docopt==0.6.2`, and `pairtools --help` succeeded after the compatibility fix.
- The OPERA-MS dedicated environment was checked on the maintainer's server: the `OPERA-MS.pl` launcher was found, `samtools --version` succeeded, `OPERA-MS.pl --help` succeeded, and the upstream `check-dependency` command completed successfully.
- `python -m unittest discover -s tests -v`: 43 tests passed.
- Python byte-compilation passed.
- CLI `--version`, top-level `--help`, and assembly-mode `run --help` passed.
- A wheel was built and inspected; required Snakefiles, environment YAMLs, and bundled HiCUP Truncater files are present, with no `.pyc` files.
- The build container provides Python 3.13, outside the declared Python 3.10–3.11 runtime range. The wheel was built and smoke-tested with `--ignore-requires-python` for packaging inspection only; installation in a supported Python 3.11 Conda environment still needs to be confirmed before release.

## End-to-end validation boundary

- A successful OPERA-MS dependency check is not the same as a complete Hybrid + Hi-C run.
- Full real-data runs of Long reads + Hi-C, Hybrid + Hi-C, and Hi-C only are not claimed by this package.
- A `--dry-run` validates the workflow plan only; it does not validate tool execution, database compatibility, or biological results.

Before production analysis, test the relevant route on the target compute system and inspect its logs, manifests, intermediate FASTA files, and final association tables.
