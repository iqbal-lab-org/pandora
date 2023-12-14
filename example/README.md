# Toy example

Here we present a walkthrough of running `pandora` on a toy example. We run:
1) `pandora` without de novo discovery;
2) `pandora` with de novo discovery (see Figure 2 of [our paper][pandora_2020_paper]).

* **There is no need to have `pandora` or `make_prg` installed. The running script will automatically download
  and run the precompiled binaries for `pandora` and `make_prg`**;

## Input data description

* `msas/` : contains the MSAs of the two genes we are using as toy example here, GC00006032 and GC00010897;
* `reads/` : contains 100x of perfect simulated reads from two toy samples. We simulated perfect reads,
where one sample genotypes to one allele of the variant sites, while the other sample genotypes towards the other allele;

## Running

```
./run_pandora.sh
```

### Quick look at the output

The output is already present in dir `out_truth`. If all went fine, the last line of the execution of the script above
should be:

```
Example run produced the expected result
```

and thus the dirs `out` and `out_truth` have the same contents.

`out/prgs`: contains output of `make_prg from_msa` and `pandora index`. Main files:
  * `pangenome.prg.fa`: the PRG itself;
  * `pangenome.prg.fa.panidx.zip`: the PRG index;
  * `pangenome.update_DS.zip`: update data structures that make the PRG updatable;

`out/pandora_discover_out`: contains the output of `pandora discover`. Main files:
  * `denovo_paths.txt`: describes the denovo paths found in all samples;

`out/updated_prgs`: contains the output of `make_prg update` and `pandora index` (on the updated PRG).
The files are similar to the ones in the `out/prgs` folder;

`out/output_toy_example_no_denovo` and `out/output_toy_example_with_denovo`: contains the output of
`pandora compare` without denovo discovery and with denovo discovery, respectively. Main files:
  * `pandora_multisample.matrix`: see https://github.com/rmcolq/pandora/wiki/FAQ#q-where-can-i-find-gene-presenceabsence-information ;
  * `pandora_multisample.vcf_ref.fa`: see https://github.com/rmcolq/pandora/wiki/FAQ#q-what-are-the-sequences-in-pandora_multisamplevcf_reffa
  * `pandora_multisample_genotyped.vcf`: the VCF file containing variants for all samples;


### Looking at the genotyped VCFs

**No denovo**

Taking a quick look at an excerpt of `out/output_toy_example_no_denovo/pandora_multisample_genotyped.vcf`
(the VCF genotyped by `pandora` without denovo sequences):

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1	toy_sample_2
GC00006032	146	.	T	C	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,78:0,96:0,78:0,96:0,157:0,192:1,0:-929.3,-67.5289:861.771	0:50,0:50,0:50,0:50,0:101,0:101,0:0,1:-17.2042,-572.517:555.313
GC00006032	160	.	A	C	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,59:0,82:0,74:0,105:0,239:0,331:1,0.25:-777.329,-53.7665:723.562	0:59,0:53,0:59,0:53,0:119,0:106,0:0,1:-24.9114,-627.779:602.868
GC00006032	218	.	T	C	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:7,28:11,44:0,28:0,44:30,57:47,88:0.75,0:-405.108,-86.4319:318.676	0:18,0:12,0:21,0:15,0:72,0:51,0:0.25,1:-23.8977,-250.155:226.257
GC00006032	237	.	T	C	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,15:0,24:0,15:0,24:0,31:0,49:1,0:-307.602,-8.43532:299.166	0:13,0:7,0:13,0:7,0:26,0:15,0:0,1:-17.8286,-204.103:186.275
GC00006032	247	.	A	C	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	.:0,0:0,0:0,0:0,0:0,0:0,0:1,1:-128,-128:0	0:10,0:4,0:10,0:4,0:10,0:4,0:0,1:-24.8363,-176.472:151.636
```

We can see samples `toy_sample_1` and `toy_sample_2` genotype towards different alleles.

**With denovo**

The VCF (`out/output_toy_example_with_denovo/pandora_multisample_genotyped.vcf`) has some new variants that were discovered and genotyped. For example:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	toy_sample_1.100x.random.illumina	toy_sample_2.100x.random.illumina
GC00006032	49	.	A	G	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	0:44,0:59,0:44,0:59,0:44,0:59,0:0,1:-15.1942,-596.333:581.138	1:0,49:0,50:0,49:0,50:0,99:0,101:1,0:-567.912,-16.6244:551.287
GC00010897	44	.	C	T	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,11:0,16:0,11:0,16:0,23:0,32:1,0:-246.34,-14.5639:231.776	0:22,0:18,0:22,0:18,0:44,0:37,0:0,1:-5.30657,-296.207:290.9
GC00010897	422	.	A	T	.	.	VC=SNP;GRAPHTYPE=SIMPLE	GT:MEAN_FWD_COVG:MEAN_REV_COVG:MED_FWD_COVG:MED_REV_COVG:SUM_FWD_COVG:SUM_REV_COVG:GAPS:LIKELIHOOD:GT_CONF	1:0,8:0,5:0,8:0,5:0,16:0,11:1,0:-181.867,-30.1108:151.756	0:12,0:9,0:12,0:9,0:12,0:9,0:0,1:-16.8478,-208.709:191.861
```

The output should be the same as using the precompiled binary.

<!--Link References-->

[pandora_2020_paper]: https://www.biorxiv.org/content/10.1101/2020.11.12.380378v2
