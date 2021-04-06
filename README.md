| Branch             | Status                                                                           |
|:-------------------|:---------------------------------------------------------------------------------|
| [`master`][master] | ![Travis (.com) branch](https://img.shields.io/travis/com/rmcolq/pandora/master) |
| [`dev`][dev]       | ![Travis (.com) branch](https://img.shields.io/travis/com/rmcolq/pandora/dev)    |

[master]: https://github.com/rmcolq/pandora/tree/master
[dev]: https://github.com/rmcolq/pandora/tree/dev


# Pandora

[TOC]: #

# Table of Contents
- [Introduction](#introduction)
- [Quick Start](#quick-start)
- [Hands-on toy example](#hands-on-toy-example)
- [Installation](#installation)
  - [Precompiled portable binary](#no-installation-needed---precompiled-portable-binary)
  - [Containers](#containers)
  - [Installation from source](#installation-from-source)
- [Usage](#usage)


## Introduction
Pandora is a tool for bacterial genome analysis using a pangenome reference graph (PanRG). It allows gene presence/absence detection and genotyping of SNPs, indels and longer variants in one or a number of samples. Pandora works with Illumina or Nanopore data. For more details, see [our paper][pandora_2020_paper].

The PanRG is a collection of 'floating'
local graphs (PRGs), each representing some orthologous region of interest
(e.g. genes, mobile elements or intergenic regions). See
https://github.com/leoisl/make_prg for a tool which can construct
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
pandora compare --genotype --illumina --max-covg 30 <panrg.fa> <read_index.tab>
```

Map Nanopore reads from a single sample to get approximate sequence for
genes present

```
pandora map <panrg.fa> <reads.fq>
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
  wget https://github.com/rmcolq/pandora/releases/download/0.9.0-rc2/pandora-linux-precompiled-v0.9.0-rc2
  ```

* **Running**:
```
chmod +x pandora-linux-precompiled-v0.9.0-rc2
./pandora-linux-precompiled-v0.9.0-rc2 -h
```

* **Notes**:
  * We provide precompiled binaries for Linux OS only;

### Containers

[![Docker Repository on Quay](https://quay.io/repository/rmcolq/pandora/status "Docker Repository on Quay")](https://quay.io/repository/rmcolq/pandora)

You can also download a containerized image of Pandora.
Pandora is hosted on Quay and images can be downloaded with the
command:

```
docker pull quay.io/rmcolq/pandora
```

Alternatively, using singularity:

```
singularity pull docker://quay.io/rmcolq/pandora
```

### Installation from source

This is the hardest way to install `pandora`, but that yields the most optimised binary.

Requirements:
- A Unix or Mac OS, with a C++11 compiler toolset (e.g. `g++`, `ld`, `make`, `ctest`, etc), `cmake`, `git` and `wget`.

- Download and install `pandora` as follows:

```
git clone --single-branch https://github.com/rmcolq/pandora.git --recursive
cd pandora
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. 
make -j4
ctest -VV
```

* If you want to produce meaningful stack traces in case `pandora` errors out, `binutils-dev` must be installed and the
  `cmake` must receive this additional parameter: `-DPRINT_STACKTRACE=True`.

## Usage

See [Usage](https://github.com/rmcolq/pandora/wiki/Usage).


<!--Link References-->
[pandora_2020_paper]: https://www.biorxiv.org/content/10.1101/2020.11.12.380378v2
