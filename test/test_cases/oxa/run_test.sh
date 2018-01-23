prg=oxa_aligned_kmeans_k15_orientated.fasta
#prg="/data2/users/rachel/projects/pandora/test/test_cases/oxa/oxa_aligned_kmeans_k15_orientated.fasta"

DATE=`date +%d_%m_%Y`

#bash ../../../header.sh &> test_index.log
#../../../build/pandora index $prg -w 5 &>> test_index.log
bash ../../../header.sh &> test_oxa9_5.15.500.$DATE.log
../../../build/pandora map -p $prg -r oxa9_read.fasta -o test_oxa9_5.15.500 -m 500 -e 0.0001 -w 5 --output_kg &>> test_oxa9_5.15.500.$DATE.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_oxa9_5.15.500*oxa*kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &> results.txt
bash ../../../header.sh &> test_JR_5.15.500.$DATE.log
../../../build/pandora map -p $prg -r JR_FAA63668_14102015_kpne_CAV1596_reads.fasta -o test_JR_5.15.500 -m 500 -w 5 --output_kg &>> test_JR_5.15.500.$DATE.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_JR_5.15.500*oxa*kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
bash ../../../header.sh &> test_all2d_JR_5.15.500.$DATE.log
perf record ../../../build/pandora map -p $prg -r JR_FAA63668_14102015_kpne_CAV1596_all2d_reads.fasta -o test_all2d_JR_5.15.500 -m 500 -w 5 --output_kg &>> test_all2d_JR_5.15.500.$DATE.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_all2d_JR_5.15.500*oxa*kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
