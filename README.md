Public Master: [![Build Status](https://travis-ci.org/rmcolq/pandora.svg?branch=master)](https://travis-ci.org/rmcolq/pandora)

Private Dev: [![Build Status](https://travis-ci.com/rmcolq/pandora.svg?token=mxzxNwUzHrkcpsL2i7zU&branch=dev)](https://travis-ci.com/rmcolq/pandora)

# pandora

### Note this is in active rapid development at present and not ready for reliable use

Pandora is a tool for bacterial genome analysis without using a reference genome,  including genetic variation from SNPs to gene presence/absence across the whole pan-genome. Core ideas are:
 - new samples look like recombinants (plus mutations) of things seen before
 - we should be analysing nucleotide-level variation everywhere, not just in core genes
 - arbitrary reference genomes are unnatural


Pandora works with Illumina or nanopore data, allowing per-sample analysis (sequence inference and SNP/indel/gene-calling) and comparison of multiple samples. To do this it uses population reference graphs (PRG) which have been built for orthologous blocks of interest (e.g. genes and intergenic regions). See https://github.com/rmcolq/make_prg for a pipeline which can construct these PRGs from a set of aligned sequence files.

It can do the following for a single sample (read dataset):

- Output inferred gene sequences for the orthologous chunks (eg genes) in the PRG
- Output a VCF showing the variation found in the pangenome genes which are present, with respect to any reference in the PRG.

For a collection of samples, it can:

- Output a matrix showing inferred copy-number of each gene in each sample genome.
- Output one VCF per orthologous-chunk, showing how samples which contained this chunk differed in their gene sequence. Variation is shown with respect to the most informative recombinant path in the PRG .
Soon, in a galaxy not so far away, it will allow

 - discovery of new variation not in the PRG


Warning - this code is still in development.

## Installation
Requires gcc 4.7 or higher on a Unix OS.

    git clone git@github.com:rmnorris/pandora.git
    cd pandora
    bash install.sh
    
## Usage
### Population Reference Graphs
Pandora assumes you have already constructed a fasta-like file of graphs, one entry for each gene/ genome region of interest. 

### Build index
Takes a fasta-like file of PRG sequences and constructs an index, and directory of gfa files to be used by pandora map.

      Usage: pandora index [options] <prgs.fa>
      Options:
      	-h,--help			Show this help message
      	-w W				Window size for (w,k)-minimizers, default 14
      	-k K				K-mer size for (w,k)-minimizers, default 15

The index stores (w,k)-minimizers for each PRG path found. These parameters can be specified, but default to w=1, k=15.

### Map reads to index
This takes a fasta of noisy long read sequence data and compares to the index. It infers which of the PRG genes/elements is present, and for those that are present it outputs the inferred sequence.

      Usage: pandora map -p PRG_FILE -r READ_FILE -o OUTDIR <option(s)>
      Options:
       -h,--help			 Show this help message
       -p,--prg_file PRG_FILE	 Specify a fasta-style prg file
       -r,--read_file READ_FILE	 Specify a file of reads in fasta format
       -o,--outdir OUTDIR	         Specify directory of output
       -w W				 Window size for (w,k)-minimizers, must be <=k, default 14
       -k K				 K-mer size for (w,k)-minimizers, default 15
       -m,--max_diff INT		 Maximum distance between consecutive hits within a cluster, default 500 (bps)
       -e,--error_rate FLOAT	 Estimated error rate for reads, default 0.11
       --genome_size NUM_BP	         Estimated length of genome, used for coverage estimation
       --output_kg			 Save kmer graphs with fwd and rev coverage annotations for found localPRGs
       --output_vcf			 Save a vcf file for each found localPRG
       --vcf_refs REF_FASTA		 A fasta file with an entry for each LocalPRG giving reference sequence for
                                     VCF. Must have a perfect match in the graph and the same name as the graph
       --illumina			 Data is from illumina rather than nanopore, so is shorter with low error rate
       --bin			 Use binomial model for kmer coverages, default is negative binomial
       --max_covg			 Maximum average coverage from reads to accept
       --regenotype			 Add extra step to carefully genotype SNP sites
      

