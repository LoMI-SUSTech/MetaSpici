# MetaSpicii
MetaSpicii (Metagenome-targeted Spatial Proximation of Inter-Chromosomal Items) is a Snakemake-based pipeline for metagenome assembly and Hi-C-guided identification of microbial genome (MAG) and mobile genetic element (MGE) associations.
MetaSpicii integrates short-read, long-read and Hi-C sequencing data into a fully automated and reproducible workflow for microbial genome reconstruction and MAG–MGE association analysis.

# Installation

## Prerequisites

Before installing MetaSpicii, please ensure that **Conda** (Miniconda or Anaconda) and **Mamba** are installed.

Mamba is recommended because it provides significantly faster dependency resolution than Conda.

If Mamba is not installed:

```bash
conda install -n base -c conda-forge mamba
```

---

## Clone MetaSpicii

Clone the repository from GitHub:

```bash
git clone <repository_url>
cd MetaSpicii
```

---

## Install MetaSpicii

Create the MetaSpicii environment using Mamba:

```bash
mamba create -n MetaSpicii \
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
conda activate MetaSpicii
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

MetaSpicii requires the following databases before running the pipeline.

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

> **Note:** Please keep track of the locations of the downloaded databases. Their paths will be required when running MetaSpicii.

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

If all commands execute successfully without errors, MetaSpicii has been installed correctly and is ready for use.

# Usage

MetaSpicii automatically detects the analysis mode based on the input files provided. No additional mode selection is required.

---

## General Usage

```bash
python MetaSpicii.py [OPTIONS]
```

Display all available options:

```bash
python MetaSpicii.py --help
```

---

## Example 1. Short-read + Hi-C

```bash
python MetaSpicii.py \
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
python MetaSpicii.py \
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
python MetaSpicii.py \
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
python MetaSpicii.py \
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
python MetaSpicii.py \
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

MetaSpicii accepts multiple sequencing strategies.

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

## Reference-guided mode

Reference-guided mode requires

- MAG directory
- vOTU fasta

Example

MAGs/
├── bin1.fa
├── bin2.fa
└── ...

vOTUs.fa
