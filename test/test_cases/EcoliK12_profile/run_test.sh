k=15
w=14
prg="pangenome_PRG.fa"

# downloard reads and extract index
wget https://s3.climb.ac.uk/nanopore/E_coli_K12_1D_R9.2_SpotON_2.pass.fasta --no-check-certificate
gunzip pangenome_PRG.fa.k15.w14.idx.gz
tar -zxvf kmer_prgs.tar.gz

# set up output dir
mkdir k$k.w$w

# run map
echo "Run pandora map"
../../../build/pandora map -p $prg -r E_coli_K12_1D_R9.2_SpotON_2.pass.fasta -o k$k.w$w/test_ecolik12 -w $w -k $k -m 250 --output_kg &>> k$k.w$w/test_ecolik12.map.log

