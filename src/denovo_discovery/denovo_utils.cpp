#include "denovo_discovery/denovo_utils.h"

std::map<std::string, std::string> get_locus_to_reads(
    const fs::path &sample_outdir,
    const std::string &sample_name
) {
    // get all reads from a locus and put in a vector (locus_to_vector_of_reads)
    std::map<std::string, std::vector<std::string>> locus_to_vector_of_reads;
    std::ifstream filtered_samfile;
    open_file_for_reading((sample_outdir / (sample_name + ".filtered.sam")).string(),
        filtered_samfile);
    std::string line;
    uint64_t segment_order = 0;
    while (std::getline(filtered_samfile, line))
    {
        std::vector<std::string> words;
        boost::split(words, line, boost::is_any_of("\t"));
        const bool is_mapped = words.size() >= 3 && (words[1] == "0" || words[1] == "16");
        if (is_mapped) {
            const std::string &read_name = words[0];
            const std::string &locus = words[2];
            const std::string &segment_seq = words[9];
            const std::string &left_flank_field = words[11];
            const std::string left_flank = left_flank_field.substr(5);
            const std::string &right_flank_field = words[12];
            const std::string right_flank = right_flank_field.substr(5);
            const std::string segment_seq_with_flanks =
                left_flank + segment_seq + right_flank;

            std::stringstream ss;
            ss << ">" << read_name << "_seg_order_" << segment_order << "\n" << segment_seq_with_flanks << "\n";
            ++segment_order;
            locus_to_vector_of_reads[locus].push_back(ss.str());
        }
    }

    // transforms the vector to a single string
    std::map<std::string, std::string> locus_to_reads;
    for(const auto &locus_to_vector_of_reads_it : locus_to_vector_of_reads) {
        const std::string &locus = locus_to_vector_of_reads_it.first;
        const std::vector<std::string> &reads = locus_to_vector_of_reads_it.second;

        std::string reads_as_a_single_string;
        for (const std::string &read : reads) {
            reads_as_a_single_string += read;
        }

        locus_to_reads[locus] = reads_as_a_single_string;
    }

    return locus_to_reads;
}


void concatenate_all_denovo_files(const std::vector<SampleData> &samples,
    const fs::path &outdir) {
    std::vector<fs::path> denovo_paths_files;
    std::vector<fs::path> denovo_sequences_files;
    for (uint32_t sample_id = 0; sample_id < samples.size(); sample_id++) {
        const auto& sample = samples[sample_id];
        const auto& sample_name = sample.first;
        fs::path denovo_path_output_file = outdir / sample_name / "denovo_paths.txt";
        denovo_paths_files.push_back(denovo_path_output_file);
        fs::path denovo_sequence_output_file = outdir / sample_name / "denovo_sequences.fa";
        denovo_sequences_files.push_back(denovo_sequence_output_file);
    }
    fs::path denovo_paths_output_file = outdir / "denovo_paths.txt";
    std::string nb_of_samples_line(std::to_string(samples.size()) + " samples");
    concatenate_text_files(denovo_paths_output_file, denovo_paths_files, nb_of_samples_line);

    fs::path denovo_sequences_output_file = outdir / "denovo_sequences.fa";
    concatenate_text_files(denovo_sequences_output_file, denovo_sequences_files);
}

void write_denovo_header_file(const std::string &sample_name,
    const fs::path &denovo_outdir, uint32_t number_of_loci_with_denovo_variants) {
    std::ofstream denovo_paths_out_header_file;
    open_file_for_writing((denovo_outdir / "denovo_paths_header.txt").string(),
        denovo_paths_out_header_file);
    denovo_paths_out_header_file << "Sample " << sample_name << std::endl;
    denovo_paths_out_header_file << number_of_loci_with_denovo_variants <<
        " loci with denovo variants" << std::endl;
    denovo_paths_out_header_file.close();
}
