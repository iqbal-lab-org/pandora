# Toy example

Here we present a walkthrough of running `pandora` on a toy example. We run: 1) `pandora` without de novo discovery; 2) [`pandora` workflow](https://github.com/iqbal-lab-org/pandora_workflow), which runs `pandora` with and without de novo discovery (see Figure 2 of [our paper][pandora_2020_paper]). Although method 2) runs both modes of `pandora`, it is much more involved than method 1), as the user needs to configure and run a `snakemake` pipeline. If no de novo discovery is required, method 1 is a lot simpler to run (two commands) as opposed to running the pipeline. For completeness, we show both methods.

## Input data description

```
msas/ : contains the MSAs of the two genes we are using as toy example here, GC00006032 and GC00010897;
prgs/toy_prg.fa : contains the PRGs of the two genes, GC00006032 and GC00010897. These are fairly simple PRGs. GC00006032 contains 4 variant sites, each representing a SNP, while GC00010897 contains 5;
reads/ : contains 100x of perfect simulated reads from two toy samples. We simulated perfect reads, where one sample genotypes to one allele of the variant sites, while the other sample genotypes towards the other allele;
pandora_workflow_data/ : contains other input and configuration files to run the pandora workflow;
```

##  pandora without de novo discovery

### Dependencies

* [`singularity`](https://sylabs.io/)

### Running
```
cd toy_example && ./run_pandora_nodenovo.sh
```

### Quick look at the output

`pandora` output will be located in directory `output_toy_example_no_denovo`.

Taking a quick look at an excerpt of the genotyped VCF (`output_toy_example_no_denovo/pandora_multisample_genotyped.vcf`):

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1	toy_sample_2
GC00006032	146	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	0:41,0:52,0:41,0:52,0:83,0:105,0:0,1:-18.7786,-526.281:507.502	1:0,15:0,15:0,15:0,15:0,31:0,31:1,0:-214.155,-3.53065:210.624
GC00006032	160	.	A	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:14,29:17,42:0,36:0,51:42,148:52,212:0.666667,0.2:-366.08,-159.943:206.137	0:12,0:8,0:18,0:9,0:38,0:24,0:0.333333,1:-20.2506,-168.103:147.853
GC00006032	218	.	T	C	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:3,11:4,14:0,11:0,14:12,23:16,28:0.75,0:-182.162,-41.9443:140.217	0:11,0:5,0:13,0:6,0:44,0:21,0:0.25,1:-19.9705,-149.683:129.712
```

We can see samples `toy_sample_1` and `toy_sample_2` genotype towards different alleles.

##  pandora workflow

### Dependencies

* [`singularity`](https://sylabs.io/)
* `git`
* `python 3.6+`

### Running
```
cd toy_example && ./run_pandora_workflow.sh
```

### Quick look at the output

`pandora` workflow output will be located at dir `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_(no|with)denovo_global_genotyping`.

Files `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_nodenovo_global_genotyping/pandora_multisample_genotyped_global.vcf` and `output_toy_example_no_denovo/pandora_multisample_genotyped.vcf` both represent running `pandora` without de novo discovery and should be very similar files, just differentiating on some header lines, and on some statistics, due to slightly different versions used. For the first file, it is the version used in the paper; for the second file, the version on the `master` branch, which is more recent, with some bugs fixed.

File `pandora_workflow/output_toy_example_workflow/illumina/100x/random/compare_withdenovo_global_genotyping/pandora_multisample_genotyped_global.vcf` is the `pandora` VCF with de novo discovery, and it has some new VCF records that were discovered and genotyped. For example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1.100x.random.illumina	toy_sample_2.100x.random.illumina
GC00006032	49	.	G	A	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,59:0,44:0,59:0,44:0,59:0,44:1,0:-570.333,-26.8805:543.452	0:50,0:48,0:50,0:48,0:100,0:97,0:0,1:-28.9415,-537.307:508.365
GC00010897	44	.	C	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,16:0,11:0,16:0,11:0,32:0,23:1,0:-220.34,-8.03511:212.304	0:18,0:22,0:18,0:22,0:37,0:44,0:0,1:-2.87264,-270.207:267.334
GC00010897	422	.	A	T	.	.	SVTYPE=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,5:0,8:0,5:0,8:0,11:0,16:1,0:-155.867,-20.2266:135.641	0:9,0:12,0:9,0:12,0:9,0:12,0:0,1:-9.39494,-182.709:173.314
```


<!--Link References-->
[pandora_2020_paper]: https://www.biorxiv.org/content/10.1101/2020.11.12.380378v2
