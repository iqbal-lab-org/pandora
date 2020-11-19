#!/usr/bin/env bash
set -eu

git clone --single-branch https://github.com/iqbal-lab-org/pandora_workflow
cp -vr msas pandora_workflow/
cp -vr prgs pandora_workflow/
cp -vr reads pandora_workflow/
cp -vr pandora_workflow_data/assemblies pandora_workflow/
cp -vr pandora_workflow_data/config pandora_workflow/
cd pandora_workflow
./setup.sh
set -eu
source venv/bin/activate
set -eu
snakemake --snakefile Snakefile_map_with_denovo        --use-singularity --configfile config/config.yaml -j1
sleep 5 # avoid locked directory issues
snakemake --snakefile Snakefile_get_denovo_updated_prg --use-singularity --configfile config/config.yaml -j1
sleep 5 # avoid locked directory issues
snakemake --snakefile Snakefile_compare                --use-singularity --configfile config/config.yaml -j1
exit 0
