# Third-Party Notices

This release includes the bundled HiCUP Truncater files and assembly-mode setup. Most scientific tools are installed as Conda packages at setup time and are not copied into the MetaSpici wheel.

## Bundled component

### HiCUP Truncater v0.9.2

- Location: `metaspici/workflow/third_party/hicup/`
- Files: `hicup_truncater`, `hicup_module.pm`, `COPYING`, `NOTICE`
- Distribution status: bundled in the MetaSpici wheel/source archive
- License: GPL-3.0-or-later, as identified in the bundled upstream notices
- Modification status: retained from the established reference workflow; runtime invocation is through Perl

## Dependencies installed or fetched at setup

The following components are installed/fetched for the assembly profile. Their source licenses continue to apply; they are not relicensed by MetaSpici.

### BASALT v1.1.0

- Source: https://github.com/PKU-EMBL/BASALT/tree/v1.1.0
- Version pin: Git tag `v1.1.0`
- Used for: short-read binning and refinement in the original short-read + Hi-C workflow
- Distribution status: source is cloned and its `BASALT_script.zip` is installed into the user-selected Conda prefix during setup; it is not embedded in the MetaSpici wheel
- License: upstream repository identifies MIT; retain the upstream license and notices with the installed source cache
- Model files: downloaded separately by the upstream BASALT model downloader into the configured cache directory; their redistribution rights should be assessed independently and they are not included in this wheel

### MetaTOR

- Source: https://github.com/koszullab/metaTOR
- Installed as a Conda package for the assembly profile
- Used for Hi-C binning, network/partition construction and MGE linking
- Upstream package metadata identifies GPL-3.0-or-later; the installed package retains its license and notices

### geNomad

- Source: https://github.com/apcamargo/genomad
- Installed as a Conda package for the assembly profile
- Used for virus/plasmid candidate detection
- Not bundled in the MetaSpici wheel; refer to the installed package's exact-version license metadata

### Other Conda-installed tools

The workflow also uses Conda-installed tools including MEGAHIT, Flye, metaMDBG, QuickMerge, seqmagick, fastp, FastUniq, dRep, CheckM, CheckV, BLAST+, BWA, samtools, samblaster, sambamba, bedtools, MetaXa2, Barrnap, tRNAscan-SE and Snakemake. They are not bundled into the wheel. Their individual package/source licenses and notices apply.

## No relicensing of third-party materials

The MetaSpici license applies only to material that the MetaSpici licensor has authority to license. Third-party tools, databases and models remain governed by their respective terms. Nothing in `LICENSE`, `COMMERCIAL_LICENSE.md`, or a MetaSpici commercial agreement overrides third-party terms unless the relevant rights holder expressly agrees in writing.

Before publishing a release, record the resolved Conda package versions, BASALT source commit, model archive identity/checksum, database versions, and the license metadata for the exact installed package set. Do not redistribute third-party databases or model files unless their licenses permit it.


## OPERA-MS (downloaded at setup time)

Hybrid + Hi-C setup downloads the official `CSB5/OPERA-MS` source repository at the `OPERA-MS-0.9.0` ref and builds it locally. OPERA-MS source and its bundled components retain their upstream license notices; users should consult the `LICENSE` and notices in the downloaded OPERA-MS tree. MetaSpici does not redistribute a prebuilt OPERA-MS binary in the wheel.
