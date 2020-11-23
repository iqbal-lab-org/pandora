#!/usr/bin/env bash
set -eu
pandora_command="pandora index prgs/toy_prg.fa && pandora compare --genotype -o output_toy_example_no_denovo prgs/toy_prg.fa reads/read_index.tsv"
singularity exec docker://rmcolq/pandora:latest bash -c "${pandora_command}"
