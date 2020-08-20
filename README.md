[![Build Status](https://travis-ci.org/rmcolq/pandora.svg?branch=master)](https://travis-ci.org/rmcolq/pandora) master

[![Build Status](https://travis-ci.com/rmcolq/pandora.svg?token=mxzxNwUzHrkcpsL2i7zU&branch=dev)](https://travis-ci.com/rmcolq/pandora) dev

![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/rmcolq/pandora)

# Pandora

## Contents
* [Introduction](#introduction)
* [Quick Start](#quick-start)
* [Installation](#installation)
* [Usage](#usage)

## Introduction
Pandora is a tool for bacterial genome analysis using a pangenome reference graph (PanRG). It allows gene presence/absence detection and genotyping of SNPs, indels and longer variants in one or a number of samples. Pandora works with Illumina or Nanopore data. Core ideas behind the method are:
 - new genomes look like recombinants (plus mutations) of things seen before
 - we should be analysing nucleotide-level variation everywhere, not just in core genes
 - arbitrary single reference genomes are unnatural and limit comparisons of diverse sets of genomes

The pangenome reference graph (PanRG) is a collection of 'floating' local graphs, each representing some orthologous region of interest (e.g. genes, mobile elements or intergenic regions). See https://github.com/rmcolq/make_prg for a pipeline which can construct these PRGs from a set of aligned sequence files.

Pandora can do the following for a single sample (read dataset):
- Output inferred mosaic of reference sequences for loci (eg genes) from the PanRG which are present
- Output a VCF showing the variation found within these loci, with respect to any reference path in the PRG.

Soon, in a galaxy not so far away, it will allow:
- discovery of new variation not in the PRG

For a collection of samples, it can:
- Output a matrix showing inferred copy-number of each locus in each sample genome
- Output a multisample pangenome VCF showing how including genotype calls for each sample in each of the loci
- Output one VCF per orthologous-chunk, showing how samples which contained this chunk differed in their gene sequence. Variation is shown with respect to the most informative recombinant path in the PRG.

Warning - this code is still in development.

## Quick Start
Index PanRG file:
```
pandora index -t 8 <panrg.fa>
```
Compare first 30X of each Illumina sample to get pangenome matrix and VCF
```
pandora compare -p <panrg.fa> -r <read_index.tab> --genotype --illumina --max_covg 30
```
Map Nanopore reads from a single sample to get approximate sequence for genes present
```
pandora map -p <panrg.fa> -r <reads.fq>
```

## Installation

### Containers
We highly recommend that you download a containerized image of Pandora. Pandora is hosted on Dockerhub and images can be downloaded with the command:
```
docker pull rmcolq/pandora:latest
```
Alternatively, using singularity:
```
singularity pull docker://rmcolq/pandora:latest
```
NB For consistency, we no longer maintain images on singularity hub.

### Installation from source
This is not recommended because the required zlib and boost system installs do not always play nicely.
If you want to take the risk:
- Requires a Unix or Mac OS.
- Requires a system install of `zlib`. If this is not already installed, [this](https://geeksww.com/tutorials/libraries/zlib/installation/installing_zlib_on_ubuntu_linux.php) tutorial is helpful or try the following.
```
wget http://www.zlib.net/zlib-1.2.11.tar.gz -O - | tar xzf -
cd zlib-1.2.11
./configure [--prefix=/prefix/path]
make
make install
```
- Requires a system installation of `boost` containing the `system`, `filesystem`, `log` (which also depends on `thread` and `date_time`) and `iostreams` libraries. If not already installed use the following or look at [this](https://www.boost.org/doc/libs/1_62_0/more/getting_started/unix-variants.html) guide.
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
### Population Reference Graphs
Pandora assumes you have already constructed a fasta-like file of graphs, one entry for each gene/ genome region of interest.
If you haven't, you will need a multiple sequence alignment for each graph. Precompiled collections of MSA representing othologous gene clusters for a number of species can be downloaded from [here](http://pangenome.de/) and converted to graphs using the pipeline from [here](https://github.com/rmcolq/make_prg).

### Build index
Takes a fasta-like file of PanRG sequences and constructs an index, and a directory of gfa files to be used by `pandora map` or `pandora compare`. These are output in the same directory as the PanRG file.

```
$ pandora index --help
Usage: pandora index [options] <prgs.fa>
Options:
        -h,--help                       Show this help message
        -w W                            Window size for (w,k)-minimizers, default 14
        -k K                            K-mer size for (w,k)-minimizers, default 15
        -t T                            Number of threads, default 1
        --offset                        Offset for PRG ids, default 0
        --outfile                       Filename for index
        --log_level                     debug,[info],warning,error
```

The index stores (w,k)-minimizers for each PanRG path found. These parameters can be specified, but default to w=14, k=15.

### Map reads to index
This takes a fasta/q of Nanopore or Illumina reads and compares to the index. It infers which of the PanRG genes/elements is present, and for those that are present it outputs the inferred sequence and a genotyped VCF.

      Usage: pandora map -p PanRG_FILE -r READ_FILE -o OUTDIR <option(s)>
      Options:
       -h,--help                        Show this help message
       -p,--prg_file PanRG_FILE         Specify a fasta-style PanRG file
       -r,--read_file READ_FILE         Specify a file of reads in fasta/q format
       -o,--outdir OUTDIR               Specify directory of output
       -w W                             Window size for (w,k)-minimizers, must be <=k, default 14
       -k K                             K-mer size for (w,k)-minimizers, default 15
       -m,--max_diff INT                Maximum distance between consecutive hits within a cluster, default 250 bps
       -e,--error_rate FLOAT            Estimated error rate for reads, default 0.11/0.001 for Nanopore/Illumina
       -c,--min_cluster_size INT        Minimum number of hits in a cluster to consider a locus present, default 10
       --genome_size NUM_BP             Estimated length of genome, used for coverage estimation, default 5000000
       --vcf_refs REF_FASTA             A fasta file with an entry for each loci in the PanRG in order, giving 
                                        reference sequence to be used as VCF ref. Must have a perfect match to a 
                                        path in the graph and the same name as the locus in the graph.
       --illumina                       Data is from Illumina, not Nanopore, so is shorter with low error rate
       --bin                            Use binomial model for kmer coverages, default is negative binomial
       --max_covg INT                   Maximum average coverage from reads to accept, default first 300
       --genotype                       Output a genotyped VCF
       --discover                       Add denovo discovery
       --denovo_kmer_size INT           Kmer size to use for denovo discovery, default 11
       --log_level LEVEL                Verbosity for logging, use "debug" for more output

### Compare reads from several samples
This takes Nanopore or Illumina read fasta/q for a number of samples, mapping each to the index. It infers which of the PanRG genes/elements is present in each sample, and outputs a presence/absence pangenome matrix, the inferred sequences for each sample and a genotyped multisample pangenome VCF.

      Usage: pandora compare -p PanRG_FILE -r READ_INDEX -o OUTDIR <option(s)>
      Options:
       -h,--help                        Show this help message
       -p,--prg_file PanRG_FILE         Specify a fasta-style PanRG file
       -r,--read_index READ_INDEX       Specify a tab delimited file with a line per sample, detailing sample id 
                                        and read fasta/q
       -o,--outdir OUTDIR               Specify directory of output
       -w W                             Window size for (w,k)-minimizers, must be <=k, default 14
       -k K                             K-mer size for (w,k)-minimizers, default 15
       -m,--max_diff INT                Maximum distance between consecutive hits within a cluster, default 250 bps
       -e,--error_rate FLOAT            Estimated error rate for reads, default 0.11/0.001 for Nanopore/Illumina
       -c,--min_cluster_size INT        Minimum number of hits in a cluster to consider a locus present, default 10
       --genome_size NUM_BP             Estimated length of genome, used for coverage estimation, default 5000000
       --vcf_refs REF_FASTA             A fasta file with an entry for each loci in the PanRG in order, giving 
                                        reference sequence to be used as VCF ref. Must have a perfect match to a 
                                        path in the graph and the same name as the locus in the graph.
       --illumina                       Data is from Illumina, not Nanopore, so is shorter with low error rate
       --bin                            Use binomial model for kmer coverages, default is negative binomial
       --max_covg INT                   Maximum average coverage from reads to accept, default first 300
       --genotype                       Output a genotyped VCF
       --log_level LEVEL                Verbosity for logging, use "debug" for more output


