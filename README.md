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

* https://github.com/PKU-EMBL/BASALT

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

## Verify Installation

After installation, verify that the required software is available:

```bash
snakemake --version
fastp --version
megahit --version
flye --version
metator --help
genomad --help
```

If the commands execute successfully, MetaSpicii has been installed correctly.

---

## Clone MetaSpicii

Clone the repository from GitHub:

```bash
git clone <repository_url>
cd MetaSpicii
```
