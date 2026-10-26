# MetaSpici V1.1 release notes

## Included

- Unified command-line entry point: `metaspici setup`, `metaspici run`, `metaspici config`, and `metaspici doctor`.
- Reference mode retained from the user-validated reference workflow baseline.
- Four assembly routes: Short reads + Hi-C, Long reads + Hi-C, Hybrid + Hi-C, and Hi-C only.
- Configurable Conda environment root and database paths.
- Dedicated OPERA-MS environment and automatic build/dependency check for Hybrid + Hi-C.
- First-time assembly setup can bootstrap the configured CheckM database path while creating a new `genome_qc` environment; an existing/incomplete environment still requires a valid database.
- Organized assembly-stage outputs, run configuration, logs, and output manifests.
- Public repository metadata, citation, license, patent notice, and third-party notices consolidated at the repository root.

## Scientific defaults retained for Reference mode

- `--min-depth 0.15`
- Strict `FE > 50`
- `--min-contacts 4`
- Optional reference-mode dRep/CheckM dereplication remains disabled by default.
- `--identity` accepts 95–99 percent.

## Validation

See `VALIDATION.md` for checks confirmed by the maintainer and the routes that still require a complete real-data pilot.
