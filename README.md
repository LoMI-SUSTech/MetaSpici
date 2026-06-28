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
