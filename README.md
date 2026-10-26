<p align="center">
<img width="220" alt="MetaSpici logo" src="https://github.com/user-attachments/assets/e031c627-7d67-4670-86ee-eecf8397328b">
</p>

# MetaSpici V1.1

**MetaSpici** (*Metagenome-targeted Spatial Proximation of
Inter-Chromosomal Items*) is a Snakemake-based pipeline for metagenome
assembly and Hi-C-guided identification of associations between
microbial genomes (MAGs) and mobile genetic elements (MGEs).

MetaSpici provides two analysis modes:

- **Assembly mode:** reconstruct MAG/MGE candidates from sequencing
  reads and infer MGE–host associations using Hi-C data.
- **Reference mode:** infer MGE–host associations from user-provided
  genomes and MGE sequences.

## Features

- Supports four assembly input configurations: Short reads + Hi-C, Long
  reads + Hi-C, Hybrid + Hi-C, and Hi-C only.
- Supports reference-guided MGE–host association analysis.
- Identifies MAG and MGE candidates from sequencing data.
- Uses Hi-C contact evidence and fold enrichment (FE) to identify
  enriched MGE–host associations.
- Provides intermediate assembly results, quality-control reports, and
  final MGE–host association tables.

## Workflow overview

| Analysis mode          | Input data                                                           | Main workflow                                                                                                               |
|------------------------|----------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|
| **Short reads + Hi-C** | Paired-end short reads and paired-end Hi-C reads                     | MEGAHIT assembly → BASALT/MetaTOR binning → combine MAG candidates → geNomad MGE identification → Hi-C association analysis |
| **Long reads + Hi-C**  | Long reads, long-read type, and paired-end Hi-C reads                | Long-read assembly → MAG/MGE candidate processing → geNomad MGE identification → Hi-C association analysis                  |
| **Hybrid + Hi-C**      | Paired-end short reads, long reads, and paired-end Hi-C reads        | OPERA-MS hybrid assembly → MAG/MGE candidate processing → geNomad MGE identification → Hi-C association analysis            |
| **Hi-C only**          | Paired-end Hi-C reads                                                | Hi-C-based assembly and binning workflow → geNomad MGE identification → Hi-C association analysis                           |
| **Reference**          | User-provided MAGs/genomes, MGE sequences, and paired-end Hi-C reads | Hi-C read processing and mapping → background-contact filtering → FE calculation → MGE–host association table               |

The assembly workflow screens contigs of at least 5 kb for geNomad-based
virus and plasmid identification.

## Installation

### Requirements

- Linux environment (for example, an HPC or Linux server).
- Conda.
- Git and network access during installation.
- Sufficient disk space for Conda environments, databases, sequencing
  reads, and intermediate results.

### 1. Download MetaSpici

```bash
git clone https://github.com/LoMI-lab/MetaSpici.git
cd MetaSpici
```

### 2. Create the MetaSpici environment

Create the launcher environment using the Conda environment file
provided with the repository:

```bash
conda env create -n MetaSpici -f metaspici/workflow/envs/launcher.yaml
conda activate MetaSpici
```

Install MetaSpici from the repository directory:

```bash
python -m pip install --no-deps .
```

Check the command-line interface:

```bash
metaspici --help
```

### Update an existing installation

If MetaSpici is already installed in a Conda launcher environment, activate that environment and install the updated source or release wheel:

```bash
conda activate MetaSpici
cd /path/to/MetaSpici
git pull
python -m pip install --no-deps --upgrade .
metaspici --version
```

If you downloaded the wheel from a GitHub Release instead, run `python -m pip install --no-deps --upgrade ./metaspici-1.1-py3-none-any.whl` from the directory containing the wheel.

## Initialize MetaSpici

### Assembly mode

```bash
metaspici setup --mode assembly \
  --env-root /path/to/metaspici/envs \
  --checkm-db /path/to/metaspici/envs/genome_qc/checkm_data \
  --checkv-db /path/to/checkv_db \
  --genomad-db /path/to/genomad_db
```

If you already have a compatible BASALT environment, add
`--basalt-prefix /absolute/path/to/BASALT`.

### Reference mode

```bash
metaspici setup --mode reference \
  --env-root /path/to/metaspici/envs \
  --checkv-db /path/to/checkv_db
```

The database and environment paths above are examples. Replace them with
the absolute paths on your system.

## Usage

Use `metaspici run --help` or `metaspici run --mode <mode> --help` to
view the available arguments.

### 1. Short reads + Hi-C

```bash
metaspici run --mode assembly \
  --short-r1 /data/sample/short_R1.fastq.gz \
  --short-r2 /data/sample/short_R2.fastq.gz \
  --hic-r1 /data/sample/HiC_R1.fastq.gz \
  --hic-r2 /data/sample/HiC_R2.fastq.gz \
  --enzyme Sau3AI \
  --sample sample1 \
  --threads 32 \
  --identity 99 \
  --outdir /results/sample1
```

### 2. Long reads + Hi-C

Supported `--long-type` values are `nanopore`, `pacbio`, `nanoporeR10`,
and `pacbioHiFi`.

```bash
metaspici run --mode assembly \
  --long-reads /data/sample/long_reads.fastq.gz \
  --long-type nanoporeR10 \
  --hic-r1 /data/sample/HiC_R1.fastq.gz \
  --hic-r2 /data/sample/HiC_R2.fastq.gz \
  --enzyme Sau3AI \
  --sample sample1_long \
  --threads 32 \
  --identity 99 \
  --outdir /results/sample1_long
```

### 3. Hybrid + Hi-C

```bash
metaspici run --mode assembly \
  --short-r1 /data/sample/short_R1.fastq.gz \
  --short-r2 /data/sample/short_R2.fastq.gz \
  --long-reads /data/sample/long_reads.fastq.gz \
  --long-type nanoporeR10 \
  --hic-r1 /data/sample/HiC_R1.fastq.gz \
  --hic-r2 /data/sample/HiC_R2.fastq.gz \
  --enzyme Sau3AI \
  --sample sample1_hybrid \
  --threads 32 \
  --identity 99 \
  --outdir /results/sample1_hybrid
```

### 4. Hi-C only

```bash
metaspici run --mode assembly \
  --hic-r1 /data/sample/HiC_R1.fastq.gz \
  --hic-r2 /data/sample/HiC_R2.fastq.gz \
  --enzyme Sau3AI \
  --sample sample1_hic \
  --threads 32 \
  --identity 99 \
  --outdir /results/sample1_hic
```

### 5. Reference mode

Provide a directory containing MAG/genome FASTA files, a FASTA file
containing the MGE sequences to test, and paired-end Hi-C reads.

```bash
metaspici run --mode reference \
  --mag-dir /data/sample/MAGs \
  --mge-fasta /data/sample/MGEs.fa \
  --hic-r1 /data/sample/HiC_R1.fastq.gz \
  --hic-r2 /data/sample/HiC_R2.fastq.gz \
  --enzyme Sau3AI \
  --sample sample1 \
  --threads 32 \
  --identity 99 \
  --outdir /results/sample1_reference
```

### Common parameters

| Parameter                  | Description                                                                             |
|----------------------------|-----------------------------------------------------------------------------------------|
| `--mode`                   | `assembly` or `reference`                                                               |
| `--short-r1`, `--short-r2` | Paired-end short reads                                                                  |
| `--long-reads`             | Long-read FASTQ file                                                                    |
| `--long-type`              | Long-read platform/type                                                                 |
| `--hic-r1`, `--hic-r2`     | Paired-end Hi-C reads                                                                   |
| `--enzyme`                 | Restriction enzyme used for Hi-C library preparation, e.g. `Sau3AI`, `MboI`, or `DpnII` |
| `--sample`                 | Sample name                                                                             |
| `--threads`                | Number of CPU threads                                                                   |
| `--identity`               | Sequence identity threshold (95–99)                                                     |
| `--outdir`                 | Output directory                                                                        |

## Input files

- **Short reads:** paired-end FASTQ files, for example
  `short_R1.fastq.gz` and `short_R2.fastq.gz`.
- **Long reads:** a FASTQ file from a supported platform (Nanopore or
  PacBio); specify the type using `--long-type`.
- **Hi-C reads:** paired-end FASTQ files, for example `HiC_R1.fastq.gz`
  and `HiC_R2.fastq.gz`.
- **Reference genomes:** a directory containing MAG/genome FASTA files,
  required for Reference mode.
- **MGE sequences:** a FASTA file containing the mobile genetic elements
  to evaluate, required for Reference mode.

## Output files

Results are written to the directory specified by `--outdir`. The exact
files depend on the selected mode and the stages completed.

Typical output groups include:

| Output | Description |
|---|---|
| `assembly_workflow/01_assembly/` | Assembly contigs and assembly-stage files |
| `assembly_workflow/02_binning/` | MAG binning and combined/dereplicated MAG candidates |
| `assembly_workflow/03_mge_detection/` | geNomad virus/plasmid identification results |
| `assembly_workflow/04_hic_network/` | Hi-C network and partition results |
| `assembly_workflow/05_mge_host_linking/` | MGE binning and MGE–host linking intermediates |
| `qc/` | Quality-control reports |
| `logs/` | Workflow and tool logs |
| `run_info/` | Run configuration and output manifests |
| `results/` | Association tables |

### Important output files

| File or directory | Description |
|---|---|
| `assembly_workflow/01_assembly/{sample}_final.contigs.fa` | Assembled contigs for the short-read assembly route |
| `assembly_workflow/02_binning/{sample}_final_binset/` | Combined/dereplicated MAG candidate set used by the downstream analysis |
| `assembly_workflow/05_mge_host_linking/metator_mge_{sample}/mges_binned.fa` | MGE sequences passed to the shared association workflow |
| `assembly_workflow/run_info/assembly_manifest.tsv` | Assembly-stage output manifest |
| `results/{sample}.associations.unfiltered_fe.tsv` | Association table before fold-enrichment filtering |
| `results/{sample}.mge_host.filtered.tsv` | Final filtered MGE–host association table |
| `run_info/config.yaml` | Effective run configuration |
| `run_info/output_manifest.tsv` | Manifest of recorded workflow outputs |

Check the files and logs in your own run directory if a particular output is missing or a stage did not complete.

## Licensing

MetaSpici uses a dual-licensing model:

- Noncommercial use is available under the [PolyForm Noncommercial
  License 1.0.0](LICENSE).
- Use outside the permitted purposes of that license requires a separate
  written [commercial license](COMMERCIAL_LICENSE.md).

PolyForm Noncommercial is a source-available license with a
noncommercial-use restriction; it is not an OSI-approved open-source
license. Read `LICENSE` and `COMMERCIAL_LICENSE.md` before use or
redistribution.

## Patents

The method implemented by MetaSpici is associated with patent rights and
a pending Hong Kong patent application. See [PATENTS.md](PATENTS.md) for
the relevant information.

## Citation

If you use MetaSpici in research, cite the software using the version or
Git commit used. Add the associated manuscript or preprint citation when
a public bibliographic record is available.

## Third-party software

MetaSpici uses third-party tools, including BASALT, MetaTOR, MEGAHIT,
long-read assemblers, dRep, geNomad, CheckM, CheckV, and OPERA-MS. Each
tool and its database remains subject to its own license and terms. See
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Support

For bug reports or feature requests, open an issue in the [MetaSpici
GitHub repository](https://github.com/LoMI-lab/MetaSpici). Include the
command used and relevant log excerpts, but do not post private
sequencing data or credentials.

For commercial licensing and technology-transfer inquiries, contact
<fanl@sustech.edu.cn>.

## Disclaimer

MetaSpici is provided under the terms of `LICENSE` and without warranty.
Results should be interpreted and validated in the context of the study.
