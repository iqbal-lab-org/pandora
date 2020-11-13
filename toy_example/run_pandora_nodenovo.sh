#!/usr/bin/env bash
set -eu

singularity pull docker://rmcolq/pandora:latest
singularity exec ./pandora-latest.simg pandora index prgs/toy_prg.fa
singularity exec ./pandora-latest.simg pandora compare --genotype -o output_toy_example_no_denovo prgs/toy_prg.fa reads/read_index.tsv
