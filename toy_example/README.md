# Toy example

Here we present a walkthrough of running `pandora` on a toy example. We run: 1) `pandora` without de novo discovery; 2) [`pandora` workflow](https://github.com/iqbal-lab-org/pandora_workflow), which runs `pandora` with and without de novo discovery (see Figure 2 of Colquhoun et al. "Nucleotide-resolution bacterial pan-genomics with reference graphs"). Although method 2) runs both modes of `pandora`, it is much more involved than method 1), as the user needs to configure and run a `snakemake` pipeline. If no de novo discovery is required, method 1 is a lot simpler to run (two commands) as opposed to running the pipeline. For completeness, we show both methods.

## Input data description

```
msas/ : contains the MSAs of the two genes we are using as toy example here, GC00006032 and GC00010897;
prgs/toy_prg.fa : contains the PRGs of the two genes, GC00006032 and GC00010897. These are fairly simple PRGs. GC00006032 contains 4 variant sites, each representing a SNP, while GC00010897 contains 5;
reads/ : contains perfect simulated reads from two toy samples. We simulated perfect reads, where one sample genotypes to one allele of the variant sites, while the other sample genotypes towards the other allele;
pandora_workflow_data : contains other input and configuration files to run the pandora workflow;
```

##  pandora without de novo discovery

### Dependencies

* `singularity`

### Running
`cd toy_example && ./run_pandora_nodenovo.sh`

### Quick look at the output

`pandora` output will be located at dir `output_toy_example_no_denovo`.

Taking a quick look at an excerpt of the genotyped VCF (`output_toy_example_no_denovo/pandora_multisample_genotyped.vcf`):

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1	toy_sample_2
GC00006032	146	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	0:6,0:11,0:6,0:11,0:12,0:22,0:0,1:-7.42551,-92.2879:84.8624	1:0,25:0,22:0,25:0,22:0,51:0,45:1,0:-312.443,-2.85628:309.587
GC00006032	160	.	A	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:2,3:3,9:0,5:0,11:6,19:11,47:0.666667,0.2:-61.987,-28.0629:33.9241	0:11,0:10,0:14,0:11,0:35,0:30,0:0.333333,1:-28.0849,-192.709:164.624
GC00006032	218	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,0:1,4:0,0:0,4:0,0:4,8:0.75,0:-28.725,-7.0005:21.7245	0:8,0:6,0:9,0:8,0:33,0:27,0:0.25,1:-30.9944,-160.472:129.478
```

We can see samples `toy_sample_1` and `toy_sample_2` genotype towards different alleles.

##  pandora workflow

### Dependencies

* `singularity`
* `git`
* `python 3.6+`

### Running
`cd toy_example && ./run_pandora_workflow.sh`

### Quick look at the output

`pandora` workflow output will be located at dir `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_(no|with)denovo_global_genotyping`.

Files `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_nodenovo_global_genotyping/pandora_multisample_genotyped_global.vcf` and `output_toy_example_no_denovo/pandora_multisample_genotyped.vcf` both represent running `pandora` without de novo discovery and should be very similar files (just differentiating on some header lines, and on some statistics, due to slightly different versions used. For the first file, it is the version used in the paper; for the second file, the version on the `master` branch, which is more recent, with some bugs fixed). 

File `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_withdenovo_global_genotyping/pandora_multisample_genotyped_global.vcf` is the `pandora` VCF with de novo discovery, and it has some new VCF records that were discovered and genotyped. For example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1.100x.random.illumina	toy_sample_2.100x.random.illumina
GC00006032	49	.	G	A	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,9:0,9:0,9:0,9:0,9:0,9:1,0:-96.8931,-8.36997:88.5231	0:69,0:58,0:69,0:58,0:139,0:116,0:0,1:-36.3338,-696.857:660.523
GC00010897	44	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,2:0,2:0,2:0,2:0,4:0,4:1,1:-32.4207,-9.39441:23.0263	0:17,0:24,0:17,0:24,0:34,0:49,0:0,1:-4.99479,-300.812:295.817
```
