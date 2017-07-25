Public Master: [![Build Status](https://travis-ci.org/rmnorris/pandora.svg?branch=master)](https://travis-ci.org/rmnorris/pandora)

Private Master: [![Build Status](https://travis-ci.com/rmnorris/pandora.svg?token=mxzxNwUzHrkcpsL2i7zU&branch=master)](https://travis-ci.com/rmnorris/pandora)

Private Dev: [![Build Status](https://travis-ci.com/rmnorris/pandora.svg?token=mxzxNwUzHrkcpsL2i7zU&branch=dev)](https://travis-ci.com/rmnorris/pandora)

# pandora
Can infer the pan-genome elements present in a sample and output highly accurate sequences for each element from long noisy reads such as Oxford Nanopore sequence data. 

Warning - this code is still in development.

## Installation
Requires gcc 4.7 or higher on a Unix OS.

    git clone pandora
    cd pandora
    bash install.sh
    
## Usage
### Build index
Takes a fasta-like file of PRG sequences and constructs an index, and directory of gfa files to be used by pandora map.

      Usage: pandora index [options] <prgs.fa>
      Options:
      	-h,--help			Show this help message
      	-w W				Window size for (w,k)-minimizers, default 1
      	-k K				K-mer size for (w,k)-minimizers, default 15

The index stores (w,k)-minimizers for each PRG path found. These parameters can be specified, but default to w=1, k=15.

### Map reads to index
This takes a fasta of noisy long read sequence data and compares to the index. It infers which of the PRG elements is present, and for those that are present it outputs the inferred sequence.

      Usage: pandora map -p PRG_FILE -r READ_FILE -o OUT_PREFIX <option(s)>
      Options:
      	-h,--help			Show this help message
      	-p,--prg_file PRG_FILE		Specify a fasta-style prg file
      	-r,--read_file READ_FILE	Specify a file of reads in fasta format
      	-o,--out_prefix OUT_PREFIX	Specify prefix of output
      	-w W				Window size for (w,k)-minimizers, default 1
      	-k K				K-mer size for (w,k)-minimizers, default 15
      	-m,--max_diff INT		Maximum distance between consecutive hits within a cluster, 
                                            default 500 (bps)
      	-e,--error_rate FLOAT		Estimated error rate for reads, default 0.11. This is adjusted 
                                            during runtime if there is sufficient coverage to infer directly from 
                                            the sequence data
      	--output_kg			Save kmer graphs with fwd and rev coverage annotations for 
                                            found localPRGs
      	--output_vcf			Save a vcf file for each found localPRG
      	--method			Method for path inference, can be max likelihood (default), 'min' 
                                            to maximize the min probability on the path, or 'both' to create outputs 
                                            with both methods
