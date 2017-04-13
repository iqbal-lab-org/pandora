profile_prefix=$1
ref=$2
fa=$3
n=$4

#make gene/ref splice
#mkdir refs
#for i in {1..266}
#    do
#	head -n $(( 2 * $i )) $3 | tail -n 2 &> refs/ref_$i.fasta
#	tail -n +2 $2 &>> refs/ref_$i.fasta
#	#awk -v ORS= '/^>/ { $0 = (NR==1 ? "" : RS) $0 RS } END { printf RS }1' refs/ref_$i.fasta &> refs/new_ref_$i.fasta
#	#mv refs/new_ref_$i.fasta refs/ref_$i.fasta 
#    done

#simulate reads
#mkdir simulated_reads
seq 33 266 | /apps/well/parallel/20161122/bin/parallel --gnu -j8 "/data2/apps/NanoSim/src/simulator.py circular -r refs/ref_{}.fasta -c $1 -o simulated_reads/read_{}.fasta -n $n"
#for i in {1..266}
#    do
#        echo $i
#	/data2/apps/NanoSim/src/simulator.py circular -r refs/ref_$i.fasta -c $1 -o simulated_reads/read_$i.fasta -n $n
#    done

