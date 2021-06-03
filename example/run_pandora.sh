#!/usr/bin/env bash
set -eu

########################################################################################################################
# argument parsing
if [[ "$#" -gt 1 || ( "$#" -eq 1  && "$1" != "conda" ) ]] ; then
    echo "Illegal parameters."
    echo "Usage: $0 or $0 conda"
    exit 1
fi
########################################################################################################################

function download_tool {
  URL=$1
  executable=$2
  wget "${URL}" -O "${executable}"
  chmod +x "${executable}"
}

# setup tools
if [ "$#" -eq 0 ] ; then
  # not conda env
  pandora_URL="https://github.com/rmcolq/pandora/releases/download/0.9.1/pandora-linux-precompiled-v0.9.1"
  pandora_executable="./pandora-linux-precompiled-v0.9.1"
  download_tool "${pandora_URL}" "${pandora_executable}"
else
  # conda env
  pandora_executable="pandora"
fi

make_prg_URL="https://github.com/leoisl/make_prg/releases/download/v0.2.0/make_prg_0.2.0"
make_prg_executable="./make_prg_0.2.0"
download_tool "${make_prg_URL}" "${make_prg_executable}"

echo "Running pandora without denovo..."
echo "Running ${make_prg_executable} from_msa"
"${make_prg_executable}" from_msa --threads 1 --input msas/ --output_prefix out/prgs/pangenome
echo "Running ${pandora_executable} index"
"${pandora_executable}" index --threads 1 out/prgs/pangenome.prg.fa
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --threads 1 --genotype -o out/output_toy_example_no_denovo out/prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running pandora without denovo - done!"

echo "Running pandora with denovo..."
echo "Running ${pandora_executable} discover"
"${pandora_executable}" discover --threads 1 --outdir out/pandora_discover_out out/prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running ${make_prg_executable} update"
"${make_prg_executable}" update --threads 1 --update_DS out/prgs/pangenome.update_DS --denovo_paths out/pandora_discover_out/denovo_paths.txt --output_prefix out/updated_prgs/pangenome_updated
echo "Running ${pandora_executable} index on updated PRGs"
"${pandora_executable}" index --threads 1 out/updated_prgs/pangenome_updated.prg.fa
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --threads 1 --genotype -o out/output_toy_example_with_denovo out/updated_prgs/pangenome_updated.prg.fa reads/read_index.tsv
echo "Running pandora with denovo - done!"

if diff -rq -I '##fileDate.*' out out_truth ; then
    echo "Example run produced the expected result"
else
    echo "ERROR: Example run DID NOT produce the expected result"
fi
