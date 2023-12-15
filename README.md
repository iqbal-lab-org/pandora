# Pandora

![master branch badge](https://github.com/rmcolq/pandora/actions/workflows/ci.yaml/badge.svg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


[TOC]: #

# Table of Contents
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Hands-on toy example](#hands-on-toy-example)
- [Installation](#installation)
  - [Precompiled portable binary](#no-installation-needed---precompiled-portable-binary)
  - [Conda](#conda)
  - [Containers](#containers)
  - [Installation from source](#installation-from-source)
- [Usage](#usage)

> Colquhoun, R.M., Hall, M.B., Lima, L. *et al.* Pandora: nucleotide-resolution bacterial pan-genomics with reference graphs. *Genome Biol* **22,** 267 (2021). https://doi.org/10.1186/s13059-021-02473-1


## Introduction
Pandora is a tool for bacterial genome analysis using a pangenome reference graph (PanRG). It allows gene presence/absence detection and genotyping of SNPs, indels and longer variants in one or a number of samples. Pandora works with Illumina or Nanopore data. For more details, see [our paper][pandora_2020_paper].

The PanRG is a collection of 'floating'
local graphs (PRGs), each representing some orthologous region of interest
(e.g. genes, mobile elements or intergenic regions). See
[make_prg][make_prg] for a tool which can construct
these PanRGs from a set of aligned sequence files.

Pandora can do the following for a single sample (read dataset):
- Output inferred mosaic of reference sequences for loci (eg genes) from the PRGs which are present in the PanRG;
- Output a VCF showing the variation found within these loci, with respect to any reference path in the PRGs;
- Discovery of new variation not in the PanRG.

For a collection of samples, it can:
- Output a matrix showing inferred presence-absence of each locus in each sample genome;
- Output a multisample pangenome VCF including genotype calls for each sample in each of the loci. Variation is shown with respect to the most informative recombinant path in the PRGs (see [our paper][pandora_2020_paper]).

> **Warning - `pandora` is not yet a production-ready tool.** 

## Quick Start

Index PanRG file:

```
pandora index -t 8 <panrg.fa>
```

Compare first 30X of each Illumina sample to get pangenome matrix and
VCF

```
pandora compare --genotype --illumina --max-covg 30 <panidx.zip> <read_index.tab>
```

Map Nanopore reads from a single sample to get approximate sequence for
genes present

```
pandora map <panidx.zip> <reads.fq>
```

## Hands-on toy example

You can test `pandora` on a toy example following [this link](example).
**There is no need to have `pandora` installed.**

## Installation

### No installation needed - precompiled portable binary

You can use `pandora` with no installation at all by simply downloading the precompiled binary, and running it.
In this binary, all libraries are linked statically.

* **Download**:
  ```
  wget https://github.com/rmcolq/pandora/releases/download/0.12.0-alpha.0/pandora-linux-precompiled-v0.12.0-alpha.0
  ```

* **Running**:
```
chmod +x pandora-linux-precompiled-v0.12.0-alpha.0
./pandora-linux-precompiled-v0.12.0-alpha.0 -h
```

* **Notes**:
  * We provide precompiled binaries for Linux OS only;

### Other installation methods

We have several other installation methods but please note that **the methods below are just available for releases,
and not for pre-releases.**
There are several improvements and changes in recent pre-releases versions that we could not still 
unit test properly, and therefore are still marked as pre-releases.

### Conda

To install `pandora` through `conda`, run **(warning: this will install version `0.9.2`)**:
```
conda install -c bioconda pandora
```

### Containers

You can also download a containerized image of Pandora.
Pandora is hosted on Quay through Biocontainers and images
can be downloaded and run with the command **(warning: this will install version `0.9.2`)**:

```
URI="quay.io/biocontainers/pandora:0.9.2--h4ac6f70_0"
docker run -it "$URI" pandora --help
```

Alternatively, you can also download and run singularity images through the galaxy project depot
**(warning: this will install version `0.9.2`)**:

```
URI="https://depot.galaxyproject.org/singularity/pandora%3A0.9.2--h4ac6f70_0"
singularity exec "$URI" pandora --help
```

### Installation from source

This is the hardest way to install `pandora`, but that yields the most optimised binary.

Requirements:
- A Unix or Mac OS, with a C++11 compiler toolset (e.g. `g++`, `ld`, `make`, `ctest`, etc), `cmake`, `git` and `wget`.

- Download and install `pandora` as follows (this example is using `4` threads, change `4` to how many threads you want):

```
git clone --single-branch https://github.com/rmcolq/pandora.git --recursive
cd pandora
mkdir -p build
cd build
cmake -DHUNTER_JOBS_NUMBER=4 -DCMAKE_BUILD_TYPE=Release ..
make -j4
ctest -VV
```

* If you want to produce meaningful stack traces in case `pandora` errors out, `binutils-dev` must be installed and the
  `cmake` must receive this additional parameter: `-DPRINT_STACKTRACE=True`.

## Usage

See [Usage](https://github.com/rmcolq/pandora/wiki/Usage).


<!--Link References-->
[pandora_2020_paper]: https://doi.org/10.1186/s13059-021-02473-1
[make_prg]: https://github.com/iqbal-lab-org/make_prg/