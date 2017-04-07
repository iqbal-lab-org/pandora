../../../build/pandora index oxa_kmeans_k15_ordered_prg.fasta -w 5 &> test_index.log
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r oxa9_read.fasta -o test_oxa9_5.15.500 -m 500 -e 0.0001 -w 5 &> test_oxa9_5.15.500.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_oxa9_5.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &> results.txt
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_reads.fasta -o test_JR_5.15.500 -m 500 -w 5 &> test_JR_5.15.500.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_JR_5.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_all2d_reads.fasta -o test_all2d_JR_5.15.500 -m 500 -w 5 &> test_all2d_JR_5.15.500.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_all2d_JR_5.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
