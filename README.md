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
- [Installation](#installation)
  - [Containers](#containers)
  - [Installation from source](#installation-from-source)
- [Usage](#usage)
  - [Population Reference Graphs](#population-reference-graphs)
  - [Build index](#build-index)
  - [Map reads to index](#map-reads-to-index)
  - [Compare reads from several samples](#compare-reads-from-several-samples)
  - [Discover novel variants](#discover-novel-variants)


## Introduction
Pandora is a tool for bacterial genome analysis using a pangenome reference graph (PanRG). It allows gene presence/absence detection and genotyping of SNPs, indels and longer variants in one or a number of samples. Pandora works with Illumina or Nanopore data. 

The PanRG is a collection of 'floating'
local graphs, each representing some orthologous region of interest
(e.g. genes, mobile elements or intergenic regions). See
https://github.com/rmcolq/make_prg for a pipeline which can construct
these PanRGs from a set of aligned sequence files.

Pandora can do the following for a single sample (read dataset):
- Output inferred mosaic of reference sequences for loci (eg genes) from the PanRG which are present;
- Output a VCF showing the variation found within these loci, with respect to any reference path in the PRG;
- Discovery of new variation not in the PanRG.

For a collection of samples, it can:
- Output a matrix showing inferred copy-number of each locus in each sample genome;
- Output a multisample pangenome VCF showing how including genotype calls for each sample in each of the loci;
- Output one VCF per locus. Variation is shown with respect to the most informative recombinant path in the PRG (see our paper).

Warning - pandora is not yet a production-ready tool yet. 

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

## Hands on a toy example

You can test `pandora` on a toy example following [this link](toy_example).
There is no need to have pandora installed, as it is run inside containers.

## Installation

### Containers

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rmcolq/pandora)

We highly recommend that you download a containerized image of Pandora.
Pandora is hosted on Dockerhub and images can be downloaded with the
command:

```
docker pull rmcolq/pandora:latest
```

Alternatively, using singularity:

```
singularity pull docker://rmcolq/pandora:latest
```

NB For consistency, we no longer maintain images on singularity hub.

### Installation from source

This is not recommended because the required zlib and boost system
installs do not always play nicely. If you want to take the risk:
- Requires a Unix or Mac OS.
- Requires a system install of `zlib`. If this is not already installed,
  [this](https://geeksww.com/tutorials/libraries/zlib/installation/installing_zlib_on_ubuntu_linux.php)
  tutorial is helpful or try the following.

```
wget http://www.zlib.net/zlib-1.2.11.tar.gz -O - | tar xzf -
cd zlib-1.2.11
./configure [--prefix=/prefix/path]
make
make install
```

- Requires a system installation of `boost` containing the `system`,
  `filesystem`, `log` (which also depends on `thread` and `date_time`)
  and `iostreams` libraries. If not already installed use the following
  or look at
  [this](https://www.boost.org/doc/libs/1_62_0/more/getting_started/unix-variants.html)
  guide.

```
wget https://sourceforge.net/projects/boost/files/boost/1.62.0/boost_1_62_0.tar.gz -O - | tar xzf -
cd boost_1_62_0
./bootstrap.sh [--prefix=/prefix/path] --with-libraries=system,filesystem,iostreams,log,thread,date_time
./b2 install
```

- Download and install `pandora` as follows:

```
git clone --single-branch https://github.com/rmcolq/pandora.git --recursive
cd pandora
mkdir -p build
cd build
cmake ..
make
ctest -VV
```

## Usage

```
$ pandora --help
Pandora: Pan-genome inference and genotyping with long noisy or short accurate reads.
Usage: pandora [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit

Subcommands:
  index                       Index population reference graph (PRG) sequences.
  map                         Quasi-map reads to an indexed PRG, infer the sequence of present loci in the sample, and (optionally) genotype/discover variants.
  compare                     Quasi-map reads from multiple samples to an indexed PRG, infer the sequence of present loci in each sample, and call variants between the samples.
  walk                        Outputs a path through the nodes in a PRG corresponding to the either an input sequence (if it exists) or the top/bottom path
  seq2path                    For each sequence, return the path through the PRG
  get_vcf_ref                 Outputs a fasta suitable for use as the VCF reference using input sequences
  random                      Outputs a fasta of random paths through the PRGs
  merge_index                 Allows multiple indices to be merged (no compatibility check
```

### Population Reference Graphs

Pandora assumes you have already constructed a fasta-like file of
graphs, one entry for each gene/ genome region of interest. If you
haven't, you will need a multiple sequence alignment for each graph.
Precompiled collections of MSA representing othologous gene clusters for
a number of species can be downloaded from [here](http://pangenome.de/)
and converted to graphs using the pipeline from
[here](https://github.com/rmcolq/make_prg).

### Build index

Takes a fasta-like file of PanRG sequences and constructs an index, and
a directory of gfa files to be used by `pandora map` or `pandora
compare`. These are output in the same directory as the PanRG file.

```
$ pandora index --help
Index population reference graph (PRG) sequences.
Usage: pandora index [OPTIONS] <PRG>

Positionals:
  <PRG> FILE [required]       PRG to index (in fasta format)

Options:
  -h,--help                   Print this help message and exit
  -w INT                      Window size for (w,k)-minimizers (must be <=k) [default: 14]
  -k INT                      K-mer size for (w,k)-minimizers [default: 15]
  -t,--threads INT            Maximum number of threads to use [default: 1]
  -o,--outfile FILE           Filename for the index [default: <PRG>.kXX.wXX.idx]
  -v                          Verbosity of logging. Repeat for increased verbosity
```

The index stores (w,k)-minimizers for each PanRG path found. These
parameters can be specified, but default to w=14, k=15.

### Map reads to index

This takes a fasta/q of Nanopore or Illumina reads and compares to the
index. It infers which of the PanRG genes/elements is present, and for
those that are present it outputs the inferred sequence and a genotyped
VCF.

```
$ pandora map --help
Quasi-map reads to an indexed PRG, infer the sequence of present loci in the sample, and optionally genotype variants.
Usage: ./pandora map [OPTIONS] <TARGET> <QUERY>

Positionals:
  <TARGET> FILE [required]    An indexed PRG file (in fasta format)
  <QUERY> FILE [required]     Fast{a,q} file containing reads to quasi-map

Options:
  -h,--help                   Print this help message and exit
  -v                          Verbosity of logging. Repeat for increased verbosity

Indexing:
  -w INT                      Window size for (w,k)-minimizers (must be <=k) [default: 14]
  -k INT                      K-mer size for (w,k)-minimizers [default: 15]

Input/Output:
  -o,--outdir DIR             Directory to write output files to [default: pandora]
  -t,--threads INT            Maximum number of threads to use [default: 1]
  --vcf-refs FILE             Fasta file with a reference sequence to use for each loci. The sequence MUST have a perfect match in <TARGET> and the same name
  --kg                        Save kmer graphs with forward and reverse coverage annotations for found loci
  --loci-vcf                  Save a VCF file for each found loci
  -C,--comparison-paths       Save a fasta file for a random selection of paths through loci
  -M,--mapped-reads           Save a fasta file for each loci containing read parts which overlapped it

Parameter Estimation:
  -e,--error-rate FLOAT       Estimated error rate for reads [default: 0.11]
  -g,--genome-size STR/INT    Estimated length of the genome - used for coverage estimation. Can pass string such as 4.4m, 100k etc. [default: 5000000]
  --bin                       Use binomial model for kmer coverages [default: negative binomial]

Mapping:
  -m,--max-diff INT           Maximum distance (bp) between consecutive hits within a cluster [default: 250]
  -c,--min-cluster-size INT   Minimum size of a cluster of hits between a read and a loci to consider the loci present [default: 10]

Preset:
  -I,--illumina               Reads are from Illumina. Alters error rate used and adjusts for shorter reads

Filtering:
  --clean                     Add a step to clean and detangle the pangraph
  --max-covg INT              Maximum coverage of reads to accept [default: 300]

Consensus/Variant Calling:
  --genotype                  Add extra step to carefully genotype sites.
  --snps                      When genotyping, only include SNP sites
  --kmer-avg INT              Maximum number of kmers to average over when selecting the maximum likelihood path [default: 100]

Genotyping:
  --local Needs: --genotype   (Intended for developers) Use coverage-oriented (local) genotyping instead of the default ML path-oriented (global) approach.
  -a INT                      Hard threshold for the minimum allele coverage allowed when genotyping [default: 0]
  -s INT                      The minimum required total coverage for a site when genotyping [default: 0]
  -D INT                      Minimum difference in coverage on a site required between the first and second maximum likelihood path [default: 0]
  -F INT                      Minimum allele coverage, as a fraction of the expected coverage, allowed when genotyping [default: 0]
  -E,--gt-error-rate FLOAT    When genotyping, assume that coverage on alternative alleles arises as a result of an error process with rate -E. [default: 0.01]
  -G,--gt-conf INT            Minimum genotype confidence (GT_CONF) required to make a call [default: 1]
```

### Compare reads from several samples

This takes Nanopore or Illumina read fasta/q for a number of samples,
mapping each to the index. It infers which of the PanRG genes/elements
is present in each sample, and outputs a presence/absence pangenome
matrix, the inferred sequences for each sample and a genotyped
multisample pangenome VCF.

```
$ pandora compare --help
Quasi-map reads from multiple samples to an indexed PRG, infer the sequence of present loci in each sample, and call variants between the samples.
Usage: ./pandora compare [OPTIONS] <TARGET> <QUERY_IDX>

Positionals:
  <TARGET> FILE [required]    An indexed PRG file (in fasta format)
  <QUERY_IDX> FILE [required] A tab-delimited file where each line is a sample identifier followed by the path to the fast{a,q} of reads for that sample

Options:
  -h,--help                   Print this help message and exit
  -v                          Verbosity of logging. Repeat for increased verbosity

Indexing:
  -w INT                      Window size for (w,k)-minimizers (must be <=k) [default: 14]
  -k INT                      K-mer size for (w,k)-minimizers [default: 15]

Input/Output:
  -o,--outdir DIR             Directory to write output files to [default: pandora]
  -t,--threads INT            Maximum number of threads to use [default: 1]
  --vcf-refs FILE             Fasta file with a reference sequence to use for each loci. The sequence MUST have a perfect match in <TARGET> and the same name
  --loci-vcf                  Save a VCF file for each found loci

Parameter Estimation:
  -e,--error-rate FLOAT       Estimated error rate for reads [default: 0.11]
  -g,--genome-size STR/INT    Estimated length of the genome - used for coverage estimation. Can pass string such as 4.4m, 100k etc. [default: 5000000]
  --bin                       Use binomial model for kmer coverages [default: negative binomial]

Mapping:
  -m,--max-diff INT           Maximum distance (bp) between consecutive hits within a cluster [default: 250]
  -c,--min-cluster-size INT   Minimum size of a cluster of hits between a read and a loci to consider the loci present [default: 10]

Preset:
  -I,--illumina               Reads are from Illumina. Alters error rate used and adjusts for shorter reads

Filtering:
  --clean                     Add a step to clean and detangle the pangraph
  --max-covg INT              Maximum coverage of reads to accept [default: 300]

Consensus/Variant Calling:
  --genotype                  Add extra step to carefully genotype sites.
  --kmer-avg INT              Maximum number of kmers to average over when selecting the maximum likelihood path [default: 100]
  
Genotyping:
  --local Needs: --genotype   (Intended for developers) Use coverage-oriented (local) genotyping instead of the default ML path-oriented (global) approach.
  -a INT                      Hard threshold for the minimum allele coverage allowed when genotyping [default: 0]
  -s INT                      The minimum required total coverage for a site when genotyping [default: 0]
  -D INT                      Minimum difference in coverage on a site required between the first and second maximum likelihood path [default: 0]
  -F INT                      Minimum allele coverage, as a fraction of the expected coverage, allowed when genotyping [default: 0]
  -E,--gt-error-rate FLOAT    When genotyping, assume that coverage on alternative alleles arises as a result of an error process with rate -E. [default: 0.01]
  -G,--gt-conf INT            Minimum genotype confidence (GT_CONF) required to make a call [default: 1]
```

### Discover novel variants

This will look for regions in the pangraph where the reads do not map
and attempt to locally assemble these regions to find novel variants.

```
$ pandora discover --help
Quasi-map reads to an indexed PRG, infer the sequence of present loci in the sample and discover novel variants.
Usage: pandora discover [OPTIONS] <TARGET> <QUERY>

Positionals:
  <TARGET> FILE [required]    An indexed PRG file (in fasta format)
  <QUERY> FILE [required]     Fast{a,q} file containing reads to quasi-map

Options:
  -h,--help                   Print this help message and exit
  --discover-k INT:[0-32)     K-mer size to use when discovering novel variants [default: 15]
  --max-ins INT               Max. insertion size for novel variants. Warning: setting too long may impair performance [default: 15]
  --covg-threshold INT        Positions with coverage less than this will be tagged for variant discovery [default: 3]
  -l INT                      Min. length of consecutive positions below coverage threshold to trigger variant discovery [default: 1]
  -L INT                      Max. length of consecutive positions below coverage threshold to trigger variant discovery [default: 30]
  -d,--merge INT              Merge candidate variant intervals within distance [default: 15]
  -N INT                      Maximum number of candidate variants allowed for a candidate region [default: 25]
  --min-dbg-dp INT            Minimum node/kmer depth in the de Bruijn graph used for discovering variants [default: 2]
  -v                          Verbosity of logging. Repeat for increased verbosity

Indexing:
  -w INT                      Window size for (w,k)-minimizers (must be <=k) [default: 14]
  -k INT                      K-mer size for (w,k)-minimizers [default: 15]

Input/Output:
  -o,--outdir DIR             Directory to write output files to [default: "pandora_discover"]
  -t,--threads INT            Maximum number of threads to use [default: 1]
  --kg                        Save kmer graphs with forward and reverse coverage annotations for found loci
  -M,--mapped-reads           Save a fasta file for each loci containing read parts which overlapped it

Parameter Estimation:
  -e,--error-rate FLOAT       Estimated error rate for reads [default: 0.11]
  -g,--genome-size STR/INT    Estimated length of the genome - used for coverage estimation. Can pass string such as 4.4m, 100k etc. [default: 5000000]
  --bin                       Use binomial model for kmer coverages [default: negative binomial]

Mapping:
  -m,--max-diff INT           Maximum distance (bp) between consecutive hits within a cluster [default: 250]
  -c,--min-cluster-size INT   Minimum size of a cluster of hits between a read and a loci to consider the loci present [default: 10]

Preset:
  -I,--illumina               Reads are from Illumina. Alters error rate used and adjusts for shorter reads

Filtering:
  --clean                     Add a step to clean and detangle the pangraph
  --clean-dbg                 Clean the local assembly de Bruijn graph
  --max-covg INT              Maximum coverage of reads to accept [default: 600]

Consensus/Variant Calling:
  --kmer-avg INT              Maximum number of kmers to average over when selecting the maximum likelihood path [default: 100]
```
