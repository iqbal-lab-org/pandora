#!/usr/bin/env bash
set -eu

# configs
pandora_URL="https://www.dropbox.com/s/3ofko0q5fcec5xe/pandora-linux-precompiled-v0.8.0-alpha?dl=1"
pandora_executable="./pandora-linux-precompiled-v0.8.0-alpha"
pandora_md5sum_file="./pandora-linux-precompiled-v0.8.0-alpha.md5sum.txt"


if md5sum -c "${pandora_md5sum_file}"; then
    # The MD5 sum match
    echo "${pandora_executable} has correct MD5 sum, proceeding..."
else
    # The MD5 sum didn't match
    echo "${pandora_executable} does not exist or does not have correct MD5 sum, downloading..."
    wget "${pandora_URL}" -O "${pandora_executable}"
    chmod +x "${pandora_executable}"
fi

"${pandora_executable}" index prgs/toy_prg.fa
"${pandora_executable}" compare --genotype -o output_toy_example_no_denovo prgs/toy_prg.fa reads/read_index.tsv

