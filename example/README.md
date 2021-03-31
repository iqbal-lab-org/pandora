# Toy example

Here we present a walkthrough of running `pandora` on a toy example. We run:
1) `pandora` without de novo discovery;
2) `pandora` with de novo discovery (see Figure 2 of [our paper][pandora_2020_paper]).

## Dependencies
* **There is no need to have `pandora` or `make_prg` installed. The running script will automatically download
  and run the precompiled binaries**;
* `wget`;
* `MAFFT` has to be in your `PATH` in order to run `make_prg update`. It can be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using conda: `conda install -c bioconda mafft`;


## Input data description

* `msas/` : contains the MSAs of the two genes we are using as toy example here, GC00006032 and GC00010897;
* `reads/` : contains 100x of perfect simulated reads from two toy samples. We simulated perfect reads, where one sample genotypes to one allele of the variant sites, while the other sample genotypes towards the other allele;

## Running

```
./run_pandora.sh
```

### Quick look at the output

`prgs`: contains output of `make_prg from_msa` and `pandora index`. Main files:
  * `pangenome.prg.fa`: the PRG itself;
  * `pangenome.prg.fa.k15.w14.idx` and `kmer_prgs`: the PRG index;
  * `pangenome.update_DS`: update data structures that make the PRG updateable;

`pandora_discover_out`: contains the output of `pandora discover`. Main files:
  * `denovo_paths.txt`: describes the denovo paths found in all samples;

`updated_prgs`: contains the output of `make_prg update` and `pandora index` (on the updated PRG).
The files are similar to the ones in the `prgs` folder;

`output_toy_example_no_denovo` and `output_toy_example_with_denovo`: contains the output of
`pandora compare` without denovo discovery and with denovo discovery, respectively. Main files:
  * `pandora_multisample.matrix`: see https://github.com/rmcolq/pandora/wiki/FAQ#q-where-can-i-find-gene-presenceabsence-information ;
  * `pandora_multisample.vcf_ref.fa`: see https://github.com/rmcolq/pandora/wiki/FAQ#q-what-are-the-sequences-in-pandora_multisamplevcf_reffa
  * `pandora_multisample_genotyped.vcf`: the VCF file containing variants for all samples;


### Looking at the genotyped VCFs

**No denovo**

Taking a quick look at an excerpt of `output_toy_example_no_denovo/pandora_multisample_genotyped.vcf` 
(the VCF genotyped by `pandora` without denovo sequences):

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1	toy_sample_2
GC00006032	146	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,41:0,52:0,41:0,52:0,83:0,105:1,0:-526.281,-18.7786:507.502	0:15,0:15,0:15,0:15,0:31,0:31,0:0,1:-3.53065,-214.155:210.624
GC00006032	160	.	A	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,26:0,40:0,33:0,50:0,106:0,160:1,0.25:-401.941,-17.9221:384.019	0:19,0:12,0:19,0:12,0:38,0:24,0:0,1:-3.32705,-218.76:215.433
GC00006032	218	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:3,11:4,14:0,11:0,14:12,23:16,28:0.75,0:-182.162,-41.9443:140.217	0:11,0:5,0:13,0:6,0:44,0:21,0:0.25,1:-19.9705,-149.683:129.712
```

We can see samples `toy_sample_1` and `toy_sample_2` genotype towards different alleles.

**With denovo**

The VCF (`output_toy_example_with_denovo/pandora_multisample_genotyped.vcf`) has some new variants that were discovered and genotyped. For example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1.100x.random.illumina	toy_sample_2.100x.random.illumina
GC00006032	49	.	A	G	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	0:44,0:59,0:44,0:59,0:44,0:59,0:0,1:-26.8805,-570.333:543.452	1:0,48:0,50:0,48:0,50:0,97:0,100:1,0:-537.307,-28.9415:508.365
GC00010897	44	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,11:0,16:0,11:0,16:0,23:0,32:1,0:-220.34,-8.03511:212.304	0:22,0:18,0:22,0:18,0:44,0:37,0:0,1:-2.87264,-270.207:267.334
GC00010897	422	.	A	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,8:0,5:0,8:0,5:0,16:0,11:1,0:-155.867,-20.2266:135.641	0:12,0:9,0:12,0:9,0:12,0:9,0:0,1:-9.39494,-182.709:173.314
```

<!--Link References-->

[pandora_2020_paper]: https://www.biorxiv.org/content/10.1101/2020.11.12.380378v2
