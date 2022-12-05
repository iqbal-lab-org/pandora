set -eu
diff <(unzip -p out/prgs/pangenome.prg.gfa.zip) <(unzip -p out_truth/prgs/pangenome.prg.gfa.zip)
diff <(unzip -p out/prgs/pangenome.update_DS.zip) <(unzip -p out_truth/prgs/pangenome.update_DS.zip)
diff <(unzip -p out/prgs/pangenome.prg.bin.zip) <(unzip -p out_truth/prgs/pangenome.prg.bin.zip)
diff <(unzip -p out/updated_prgs/pangenome_updated.prg.bin.zip) <(unzip -p out_truth/updated_prgs/pangenome_updated.prg.bin.zip)
diff <(unzip -p out/updated_prgs/pangenome_updated.update_DS.zip) <(unzip -p out_truth/updated_prgs/pangenome_updated.update_DS.zip)
diff <(unzip -p out/updated_prgs/pangenome_updated.prg.gfa.zip) <(unzip -p out_truth/updated_prgs/pangenome_updated.prg.gfa.zip)
