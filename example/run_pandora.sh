#!/usr/bin/env bash
set -eu

# configs
pandora_URL="https://github.com/rmcolq/pandora/releases/download/pandora_paper_tag1/pandora-linux-precompiled-pandora_paper_tag1"
pandora_executable="./pandora-linux-precompiled-pandora_paper_tag1"
make_prg_URL="https://github.com/leoisl/make_prg/releases/download/v0.2.0_prototype/make_prg_0.2.0_prototype"
make_prg_executable="./make_prg_0.2.0_prototype"


function download_tool {
  URL=$1
  executable=$2
  wget "${URL}" -O "${executable}"
  chmod +x "${executable}"
}

download_tool "${pandora_URL}" "${pandora_executable}"
download_tool "${make_prg_URL}" "${make_prg_executable}"

echo "Running pandora without denovo..."
echo "Running ${make_prg_executable} from_msa"
"${make_prg_executable}" from_msa --input msas/ --output_prefix prgs/pangenome
echo "Running ${pandora_executable} index"
"${pandora_executable}" index prgs/pangenome.prg.fa
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --genotype -o output_toy_example_no_denovo prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running pandora without denovo - done!"

echo "Running pandora with denovo..."
echo "Running ${pandora_executable} discover"
"${pandora_executable}" discover --outdir pandora_discover_out prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running ${make_prg_executable} update"
"${make_prg_executable}" update --update_DS prgs/pangenome.update_DS --denovo_paths pandora_discover_out/denovo_paths.txt --output_prefix updated_prgs/pangenome_updated
echo "Running ${pandora_executable} index on updated PRGs"
"${pandora_executable}" index updated_prgs/pangenome_updated.prg.fa
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --genotype -o output_toy_example_with_denovo updated_prgs/pangenome_updated.prg.fa reads/read_index.tsv
echo "Running pandora with denovo - done!"
