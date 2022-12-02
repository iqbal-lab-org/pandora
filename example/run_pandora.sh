#!/usr/bin/env bash
set -eu

########################################################################################################################
# configs
pandora_version="0.10.0-alpha.0"
pandora_URL="https://github.com/rmcolq/pandora/releases/download/${pandora_version}/pandora_${pandora_version}"
make_prg_version="0.4.0"
make_prg_URL="https://github.com/iqbal-lab-org/make_prg/releases/download/${make_prg_version}/make_prg_${make_prg_version}"
########################################################################################################################

########################################################################################################################
# argument parsing
if [[ "$#" -gt 1 || ( "$#" -eq 1  && "$1" != "conda" ) ]] ; then
    echo "Illegal parameters."
    echo "Usage: $0 or $0 conda"
    exit 1
fi
########################################################################################################################

########################################################################################################################
# setup tools
function download_tool {
  URL=$1
  executable=$2
  wget "${URL}" -O "${executable}"
  chmod +x "${executable}"
}

if [ "$#" -eq 0 ] ; then
  # not conda env
  pandora_executable="./pandora-linux-precompiled-v${pandora_version}"
  if [ ! -f ${pandora_executable} ]; then
    echo "pandora executable not found, downloading it..."
    download_tool "${pandora_URL}" "${pandora_executable}"
  fi

  make_prg_executable="./make_prg-v${make_prg_version}"
  if [ ! -f ${make_prg_executable} ]; then
    echo "make_prg executable not found, downloading it..."
    download_tool "${make_prg_URL}" "${make_prg_executable}"
  fi
else
  # conda env
  pandora_executable="pandora"
  make_prg_executable="make_prg"
fi
########################################################################################################################

########################################################################################################################
# run sample example
echo "Running pandora without denovo..."
echo "Running ${make_prg_executable} from_msa"
${make_prg_executable} from_msa --threads 1 --input msas/ --output-prefix out/prgs/pangenome
echo "Running ${pandora_executable} index"
"${pandora_executable}" index --threads 1 out/prgs/pangenome.prg.fa
echo "Running ${pandora_executable} map"
"${pandora_executable}" map --threads 1 --genotype -o out/map_toy_sample_1 out/prgs/pangenome.prg.fa reads/toy_sample_1/toy_sample_1.100x.random.illumina.fastq
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --threads 1 --genotype -o out/output_toy_example_no_denovo out/prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running pandora without denovo - done!"

echo "Running pandora with denovo..."
echo "Running ${pandora_executable} discover"
"${pandora_executable}" discover --threads 1 --outdir out/pandora_discover_out out/prgs/pangenome.prg.fa reads/read_index.tsv
echo "Running ${make_prg_executable} update"
${make_prg_executable} update --threads 1 --update-DS out/prgs/pangenome.update_DS.zip --denovo-paths out/pandora_discover_out/denovo_paths.txt --output-prefix out/updated_prgs/pangenome_updated
echo "Running ${pandora_executable} index on updated PRGs"
"${pandora_executable}" index --threads 1 out/updated_prgs/pangenome_updated.prg.fa
echo "Running ${pandora_executable} compare"
"${pandora_executable}" compare --threads 1 --genotype -o out/output_toy_example_with_denovo out/updated_prgs/pangenome_updated.prg.fa reads/read_index.tsv
echo "Running pandora with denovo - done!"

# first compare non-zip files
if diff -rq -I '##fileDate.*' --exclude "*.zip" out out_truth ; then
  # now compare zip files
  if scripts/diff_zip_files.sh ; then
      echo "Example run produced the expected result"
  else
      echo "ERROR: Example run DID NOT produce the expected result"
  fi
else
    echo "ERROR: Example run DID NOT produce the expected result"
fi
########################################################################################################################
