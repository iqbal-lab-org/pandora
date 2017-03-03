../../../build/pandora index oxa_kmeans_k15_ordered_prg.fasta &> test_index.log
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r oxa9_read.fasta -o test_oxa9_1.15.500 -m 500 -e 0.0001 &> test_oxa9_1.15.500.log
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_reads.fasta -o test_JR_1.15.500 -m 500 &> test_JR_1.15.500.log
../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_all2d_reads.fasta -o test_all2d_JR_1.15.500 -m 500 &> test_all2d_JR_1.15.500.log
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_oxa9_1.15.500_oxa_mlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &> results.txt
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_JR_1.15.500_oxa_mlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_all2d_JR_1.15.500_oxa_mlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> results.txt
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_oxa9_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &> new_results.txt
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_JR_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> new_results.txt
echo -e "$(head -2 oxa9_read.fasta | tail -n 1)\n$(head -2 test_all2d_JR_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> new_results.txt
echo -e "$(head -2 test_oxa9_1.15.500_oxa_mlp.fasta | tail -n 1)\n$(head -2 test_oxa9_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &> diff_results.txt
echo -e "$(head -2 test_JR_1.15.500_oxa_mlp.fasta | tail -n 1)\n$(head -2 test_JR_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> diff_results.txt
echo -e "$(head -2 test_all2d_JR_1.15.500_oxa_mlp.fasta | tail -n 1)\n$(head -2 test_all2d_JR_1.15.500_oxa_kmlp.fasta | tail -n 1)" | ~/apps/cortex/scripts/analyse_variants/needleman_wunsch/needleman_wunsch --stdin --zam &>> diff_results.txt
#../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r oxa9_read.fasta -o test_oxa9_1.15.500 -m 500 -e 0.0001 --output_p_dist &> test_oxa9_1.15.500.pdist.log
#../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_reads.fasta -o test_JR_1.15.500 -m 500 --output_p_dist &> test_JR_1.15.500.pdist.log
#../../../build/pandora map -p oxa_kmeans_k15_ordered_prg.fasta -r JR_FAA63668_14102015_kpne_CAV1596_all2d_reads.fasta -o test_all2d_JR_1.15.500 -m 500 --output_p_dist &> test_all2d_JR_1.15.500.pdist.log
