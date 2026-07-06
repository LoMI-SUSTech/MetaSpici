# MetaSpici
MetaSpici (Metagenome-targeted Spatial Proximation of Inter-Chromosomal Items) is a Snakemake-based pipeline for metagenome assembly and Hi-C-guided identification of microbial genome (MAG) and mobile genetic element (MGE) associations.
MetaSpici integrates short-read, long-read and Hi-C sequencing data into a fully automated and reproducible workflow for microbial genome reconstruction and MAG–MGE association analysis.

# Installation

## Prerequisites

Before installing MetaSpici, please ensure that **Conda** (Miniconda or Anaconda) and **Mamba** are installed.

Mamba is recommended because it provides significantly faster dependency resolution than Conda.

If Mamba is not installed:

```bash
conda install -n base -c conda-forge mamba
```

---

## Clone MetaSpici

Clone the repository from GitHub:

```bash
git clone <repository_url>
cd MetaSpici
```

---

## Install MetaSpici

Create the MetaSpici environment using Mamba:

```bash
mamba create -n MetaSpici \
    -c conda-forge \
    -c bioconda \
    python=3.9 \
    snakemake \
    fastp fastuniq \
    hicstuff \
    megahit flye metamdbg \
    metator drep \
    genomad \
    bwa samtools samblaster bamm \
    seqkit seqmagick \
    pandas numpy pyyaml openpyxl \
    trnascan-se \
    metaxa \
    barrnap
```

Activate the environment:

```bash
conda activate MetaSpici
```

---

## Install BASALT

Due to its complex dependency structure, BASALT is recommended to be installed in a separate Conda environment.

Please follow the official installation instructions provided by the BASALT project:

https://github.com/PKU-EMBL/BASALT

A typical installation is:

```bash
git clone https://github.com/PKU-EMBL/BASALT.git
cd BASALT

mamba env create -n BASALT -f basalt_env.yml
```

Activate the BASALT environment when required:

```bash
conda activate BASALT
```

---

## Download Required Databases

MetaSpici requires the following databases before running the pipeline.

### geNomad database

Download the geNomad database:

```bash
genomad download-database database/genomad_db
```

Alternatively, specify another directory:

```bash
genomad download-database /path/to/genomad_db
```

---

### CheckV database

Download the CheckV database:

```bash
checkv download_database database/checkv_db
```

Or specify another directory:

```bash
checkv download_database /path/to/checkv_db
```

> **Note:** Please keep track of the locations of the downloaded databases. Their paths will be required when running MetaSpici.

---

## Verify Installation

Verify that the required software has been installed successfully:

```bash
snakemake --version
fastp --version
hicstuff --help
megahit --version
flye --version
metamdbg --help
metator --help
dRep --help
genomad --help
bwa
samtools --version
seqkit version
trnascan-se -h
metaxa -h
barrnap --help
```

If all commands execute successfully without errors, MetaSpici has been installed correctly and is ready for use.

# Usage

MetaSpici automatically detects the analysis mode based on the input files provided. No additional mode selection is required.

---

## General Usage

```bash
python MetaSpici.py [OPTIONS]
```

Display all available options:

```bash
python MetaSpici.py --help
```

---

## Example 1. Short-read + Hi-C

```bash
python MetaSpici.py \
    -s1 short_R1.fastq.gz \
    -s2 short_R2.fastq.gz \
    -h1 hic_R1.fastq.gz \
    -h2 hic_R2.fastq.gz \
    -e DpnII \
    -n fecal \
    -j 32 \
    -o results/
```

---

## Example 2. Long-read + Hi-C (Oxford Nanopore)

```bash
python MetaSpici.py \
    --nanopore nanopore.fastq.gz \
    -h1 hic_R1.fastq.gz \
    -h2 hic_R2.fastq.gz \
    -e DpnII \
    -n fecal \
    -j 32 \
    -o results/
```

---

## Example 3. Hybrid Assembly (Short + Long + Hi-C)

```bash
python MetaSpici.py \
    -s1 short_R1.fastq.gz \
    -s2 short_R2.fastq.gz \
    --nanopore nanopore.fastq.gz \
    -h1 hic_R1.fastq.gz \
    -h2 hic_R2.fastq.gz \
    -e DpnII \
    -n fecal \
    -j 32 \
    -o results/
```

---

## Example 4. Hi-C Only

```bash
python MetaSpici.py \
    -h1 hic_R1.fastq.gz \
    -h2 hic_R2.fastq.gz \
    -e DpnII \
    -n fecal \
    -j 32 \
    -o results/
```

---

## Example 5. Reference-guided Analysis

```bash
python MetaSpici.py \
    --reference-guided \
    --mag_dir MAGs/ \
    --votu vOTUs.fa \
    -h1 hic_R1.fastq.gz \
    -h2 hic_R2.fastq.gz \
    -e DpnII \
    -n fecal \
    -j 32 \
    -o results/
```

---

## Command-line Arguments

| Argument             | Description                                            |
| -------------------- | ------------------------------------------------------ |
| `-s1`, `-s2`         | Short-read FASTQ files (paired-end)                    |
| `--nanopore`         | Oxford Nanopore reads                                  |
| `--nanopore_R10`     | Oxford Nanopore R10 reads                              |
| `--pacbio`           | PacBio CLR reads                                       |
| `--pacbioHiFi`       | PacBio HiFi reads                                      |
| `-h1`, `-h2`         | Hi-C paired-end FASTQ files (required)                 |
| `-e`                 | Restriction enzyme used for Hi-C library preparation   |
| `-n`                 | Sample name                                            |
| `-j`                 | Number of CPU threads                                  |
| `-o`                 | Output directory                                       |
| `--reference-guided` | Enable reference-guided MAG–vOTU association mode      |
| `--mag_dir`          | Directory containing MAG genomes                       |
| `--votu`             | vOTU FASTA file                                        |
| `--target`           | Execute a specific Snakemake target                    |
| `--create-envs-only` | Create Conda environments without running the workflow |


# Input Files

MetaSpici accepts multiple sequencing strategies.

## Short reads

Paired-end Illumina FASTQ files.

Example:

short_R1.fastq.gz
short_R2.fastq.gz

---

## Long reads

Supported platforms

- Oxford Nanopore
- Oxford Nanopore R10
- PacBio CLR
- PacBio HiFi

Example

nanopore.fastq.gz

---

## Hi-C reads

Paired-end Hi-C FASTQ files.

Example

hic_R1.fastq.gz
hic_R2.fastq.gz

---

### Reference-guided mode

The reference-guided mode is designed for users who already have **reference MAGs** and **viral genomes (vOTUs)** available. Instead of performing de novo assembly and genome binning, MetaSpici directly utilizes these reference sequences together with Hi-C sequencing data to identify MAG–MGE associations.

This mode is particularly useful when:

- High-quality MAGs have already been reconstructed from the same metagenomic sample.
- Reference genomes are available from isolates or pure-culture sequencing obtained from the same environment.
- Existing genome catalogs have been generated in previous studies and users only wish to perform Hi-C-based MAG–MGE association analysis.

**Note:** MetaSpici does **not** generate MAGs or vOTUs in reference-guided mode. These reference genomes must be prepared independently by the user before running the pipeline.

#### Required inputs

**MAG directory**

A directory containing assembled MAG genomes in FASTA format.

Example:

```text
MAGs/
├── MAG_001.fa
├── MAG_002.fa
├── MAG_003.fa
└── ...
```

**vOTU FASTA**

A FASTA file containing viral genomes or viral operational taxonomic units (vOTUs).

Example:

```text
vOTUs.fa
```

**Hi-C reads**

Paired-end Hi-C sequencing reads generated from the same sample.

```text
hic_R1.fastq.gz
hic_R2.fastq.gz
```

# Output Files

MetaSpici stores all analysis results in the user-specified output directory. A typical output structure is shown below:

```text
results/
├── assembly/
├── mapping/
├── binning/
├── blast_results_mock/
├── filtered_beds_clean/
├── hic_results/
├── mapped_reads/
├── output/
├── repeat_regions/
├── tmp/
├── merged_links.tsv
└── config.yaml
```

---

## Directory Description

| Directory/File         | Description                                                                                                                                           |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `assembly/`            | Genome assembly results, including assembled contigs generated from short reads, long reads or hybrid assembly.                                       |
| `mapping/`             | Read preprocessing and mapping results, including quality control, duplicate removal, Hi-C restriction site processing, and reference genome masking. |
| `binning/`             | Genome binning and refinement results for metagenome-assembled genomes (MAGs).                                                                        |
| `blast_results_mock/`  | BLAST self-alignment results used to identify homologous or repetitive genomic regions for downstream masking.                                        |
| `repeat_regions/`      | Genomic regions identified as repetitive or highly similar based on BLAST analysis.                                                                   |
| `filtered_beds_clean/` | Filtered Hi-C contact files after masking repetitive regions and removing contacts originating from these regions.                                    |
| `mapped_reads/`        | Intermediate files generated during read alignment and mapping.                                                                                       |
| `hic_results/`         | Final filtered Hi-C interaction data and abundance estimation results.                                                                                |
| `output/`              | Intermediate and final results generated during mobile genetic element (MGE) identification and MAG–MGE association analysis.                         |
| `tmp/`                 | Temporary files generated during pipeline execution.                                                                                                  |
| `merged_links.tsv`     | Unfiltered MAG–MGE association results generated from Hi-C interaction data.                                                                          |
| `config.yaml`          | Configuration file automatically generated from the command-line arguments.                                                                           |

---

## Primary Outputs

The major results produced by MetaSpici include:

* High-quality assembled contigs.
* Metagenome-assembled genomes (MAGs).
* Filtered Hi-C interaction datasets.
* Hi-C-based microbial abundance estimation.
* Predicted mobile genetic elements (MGEs).
* Hi-C-supported MAG–MGE association results.
