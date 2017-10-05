k=$1
w=$2
prg=oxa_aligned_kmeans_k15_orientated.fasta

mkdir k$k.w$w
#bash ../../../header.sh &> test_index.k$k.w$w.log
../../../build/pandora index $prg -k $k -w $w &>> test_index.k$k.w$w.log

seq 2 2 532 | /apps/well/parallel/20161122/bin/parallel --gnu -j8 "/data2/users/rachel/projects/pandora_dev/build/pandora map -p $prg -r reads/oxa_read_{}.fasta -o k$k.w$w/{} -w $w -k $k -m 500 -e 0.0001 &>> k$k.w$w/{}.log"

for i in {2..532..2}
    do
	echo -e "$(head -1 reads/oxa_read_$i.fasta)" &>> k$k.w$w/results.txt
        echo -e "$(head -2 reads/oxa_read_$i.fasta | tail -n 1)\n$(head -2 k$k.w$w/$i.oxa.kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> k$k.w$w/results.txt
    done
