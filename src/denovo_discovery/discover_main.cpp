//
// Created by michael on 26/9/20.
//

#include "denovo_discovery/discover_main.h"
#include <boost/algorithm/string.hpp>
#include "utils.h"
#include "denovo_discovery/racon.h"

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
    while (std::getline(filtered_samfile, line))
    {
        std::vector<std::string> words;
        boost::split(words, line, boost::is_any_of("\t"));
        const bool is_mapped = words.size() >= 3 && words[1] == "0";
        if (is_mapped) {
            const std::string &read_name = words[0];
            const std::string &locus = words[2];
            const std::string &read_seq = words[9];
            std::stringstream ss;
            ss << ">" << read_name << "\n" << read_seq << "\n";
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

void setup_discover_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<DiscoverOptions>();
    auto* discover_subcmd = app.add_subcommand("discover",
        "Quasi-map reads to an indexed PRG, infer the "
        "sequence of present loci in the sample and discover novel variants.");

    discover_subcmd
        ->add_option("<TARGET>", opt->prgfile, "An indexed PRG file (in fasta format)")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    std::string description
        = "A tab-delimited file where each line is a sample identifier followed by "
          "the path to the fast{a,q} of reads for that sample";
    discover_subcmd->add_option("<QUERY_IDX>", opt->reads_idx_file, description)
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    discover_subcmd
        ->add_option(
            "-w", opt->window_size, "Window size for (w,k)-minimizers (must be <=k)")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    discover_subcmd->add_option("-k", opt->kmer_size, "K-mer size for (w,k)-minimizers")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    discover_subcmd
        ->add_option("-o,--outdir", opt->outdir, "Directory to write output files to")
        ->transform(make_absolute)
        ->type_name("DIR")
        ->capture_default_str()
        ->group("Input/Output");

    discover_subcmd
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Input/Output");

    discover_subcmd
        ->add_option(
            "-e,--error-rate", opt->error_rate, "Estimated error rate for reads")
        ->capture_default_str()
        ->group("Parameter Estimation");

    discover_subcmd
        ->add_option("-g,--genome-size", opt->genome_size,
            "Estimated length of the genome - used for coverage estimation. Can pass "
            "string such as 4.4m, 100k etc.")
        ->transform(transform_cli_gsize)
        ->capture_default_str()
        ->type_name("STR/INT")
        ->group("Parameter Estimation");

    discover_subcmd
        ->add_option("-m,--max-diff", opt->max_diff,
            "Maximum distance (bp) between consecutive hits within a cluster")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    discover_subcmd
        ->add_flag("--kg", opt->output_kg,
            "Save kmer graphs with forward and reverse coverage annotations for found "
            "loci")
        ->group("Input/Output");

    discover_subcmd
        ->add_flag("-M,--mapped-reads", opt->output_mapped_read_fa,
            "Save a fasta file for each loci containing read parts which overlapped it")
        ->group("Input/Output");

    discover_subcmd
        ->add_flag("-I,--illumina", opt->illumina,
            "Reads are from Illumina. Alters error rate used and adjusts for shorter "
            "reads")
        ->group("Preset");

    discover_subcmd
        ->add_flag(
            "--clean", opt->clean, "Add a step to clean and detangle the pangraph")
        ->group("Filtering");

    discover_subcmd
        ->add_flag(
            "--clean-dbg", opt->clean_dbg, "Clean the local assembly de Bruijn graph")
        ->group("Filtering");

    discover_subcmd
        ->add_flag("--bin", opt->binomial,
            "Use binomial model for kmer coverages [default: negative binomial]")
        ->group("Parameter Estimation");

    discover_subcmd
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Filtering");

    discover_subcmd
        ->add_option("--discover-k", opt->denovo_kmer_size,
            "K-mer size to use when discovering novel variants")
        ->capture_default_str()
        ->check(CLI::Range(0, MAX_DENOVO_K)
                    .description("[0-" + std::to_string(MAX_DENOVO_K) + ")"))
        ->type_name("INT");

    description = "Max. insertion size for novel variants. Warning: "
                  "setting too long may impair performance";
    discover_subcmd->add_option("--max-ins", opt->max_insertion_size, description)
        ->capture_default_str()
        ->type_name("INT");

    description
        = "Positions with coverage less than this will be tagged for variant discovery";
    discover_subcmd
        ->add_option("--covg-threshold", opt->min_candidate_covg, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Min. length of consecutive positions below coverage threshold to "
                  "trigger variant discovery";
    discover_subcmd->add_option("-l", opt->min_candidate_len, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Max. length of consecutive positions below coverage threshold to "
                  "trigger variant discovery";
    discover_subcmd->add_option("-L", opt->max_candidate_len, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Merge candidate variant intervals within distance";
    discover_subcmd->add_option("-d,--merge", opt->merge_dist, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Maximum number of candidate variants allowed for a candidate region";
    discover_subcmd->add_option("-N", opt->max_num_candidate_paths, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Minimum node/kmer depth in the de Bruijn graph used for discovering "
                  "variants";
    discover_subcmd
        ->add_option(
            "--min-dbg-dp", opt->min_covg_for_node_in_assembly_graph, description)
        ->capture_default_str()
        ->type_name("INT");

    description
        = "Minimum size of a cluster of hits between a read and a loci to consider "
          "the loci present";
    discover_subcmd
        ->add_option("-c,--min-cluster-size", opt->min_cluster_size, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    description = "Maximum number of kmers to average over when selecting the maximum "
                  "likelihood path";
    discover_subcmd->add_option("--kmer-avg", opt->max_num_kmers_to_avg, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    discover_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    discover_subcmd->callback([opt]() { pandora_discover(*opt); });
}

void concatenate_denovo_files(
    const fs::path& output_filename, const std::vector<fs::path>& input_filenames)
{
    ofstream output_filehandler;
    open_file_for_writing(output_filename.string(), output_filehandler);

    output_filehandler << input_filenames.size() << " samples" << std::endl;

    for (const fs::path& input_filename : input_filenames) {
        ifstream input_filehandler;
        open_file_for_reading(input_filename.string(), input_filehandler);
        output_filehandler << input_filehandler.rdbuf();
        input_filehandler.close();
    }

    output_filehandler.close();
}

void find_denovo_variants_core(std::vector<CandidateRegion*>& candidate_regions,
    const SampleIdText& sample_name, const fs::path& sample_outdir,
    const DenovoDiscovery& denovo, uint32_t child_id, uint32_t threads)
{
    // create temp dir
    const auto temp_dir { sample_outdir / ("temp_child_" + int_to_string(child_id)) };
    fs::create_directories(temp_dir);

    CandidateRegionWriteBuffer buffer(sample_name);
    for (uint32_t candidate_region_index = child_id;
         candidate_region_index < candidate_regions.size();
         candidate_region_index += threads) {
        CandidateRegion& candidate_region { *(
            candidate_regions[candidate_region_index]) };
        denovo.find_paths_through_candidate_region(candidate_region, temp_dir);
        candidate_region.write_denovo_paths_to_buffer(buffer);
    }

    // serialise the buffer
    {
        const auto buffer_binary_filename
            = temp_dir / "candidate_regions_write_buffer.bin";
        std::ofstream buffer_binary_filehandler(buffer_binary_filename.string());
        boost::archive::text_oarchive buffer_binary_archive(buffer_binary_filehandler);
        // write class instance to archive
        buffer_binary_archive << buffer;
        // archive and stream closed when destructors are called
    }
}

void find_denovo_variants_multiprocess(CandidateRegions& candidate_regions,
    const SampleIdText& sample_name, const fs::path& sample_outdir,
    const DenovoDiscovery& denovo, uint32_t threads)
{
    // transforms CandidateRegions into a vector of pointers to CandidateRegion so that
    // it is easier to multithread/multiprocess on it
    std::vector<CandidateRegion*> candidate_regions_as_vector;
    candidate_regions_as_vector.reserve(candidate_regions.size());
    for (auto& element : candidate_regions) {
        CandidateRegion* candidate_region_pointer = &(element.second);
        candidate_regions_as_vector.push_back(candidate_region_pointer);
    }

    // forking due to GATB
    size_t child_id;
    bool on_child;
    for (child_id = 0; child_id < threads; ++child_id) {
        int child_process_id = fork();
        bool error_creating_child_process = child_process_id == -1;
        on_child = child_process_id == 0;

        if (error_creating_child_process) {
            fatal_error("Error creating child process.");
        } else if (on_child) {
            break;
        } else {
            BOOST_LOG_TRIVIAL(info) << "Child process id " << child_process_id
                                    << " (child #" << child_id << ") created...";
        }
    }

    if (on_child) {
        find_denovo_variants_core(candidate_regions_as_vector, sample_name,
            sample_outdir, denovo, child_id, threads);
        std::exit(0);
    } else {
        // wait for all children to finish
        for (child_id = 0; child_id < threads; ++child_id) {
            int child_pid = wait(NULL);
            bool error_on_waiting_for_child = child_pid == -1;
            if (error_on_waiting_for_child) {
                fatal_error("Error waiting for child process.");
            } else {
                BOOST_LOG_TRIVIAL(info)
                    << "Child process " << child_pid << " finished!";
            }
        }
    }

    // add all candidate region write buffers to a central one
    CandidateRegionWriteBuffer buffer(sample_name);
    std::vector<fs::path> buffer_binary_filenames;
    for (child_id = 0; child_id < threads; child_id++) {
        CandidateRegionWriteBuffer child_buffer;
        {
            const auto child_temp_dir { sample_outdir
                / ("temp_child_" + int_to_string(child_id)) };
            const auto child_buffer_binary_filename
                = child_temp_dir / "candidate_regions_write_buffer.bin";
            // create and open an archive for input
            std::ifstream child_buffer_binary_filehandler(
                child_buffer_binary_filename.string());
            boost::archive::text_iarchive child_buffer_binary_archive(
                child_buffer_binary_filehandler);
            // read class state from archive
            child_buffer_binary_archive >> child_buffer;
            // archive and stream closed when destructors are called
        }
        buffer.merge(child_buffer);
    }

    auto denovo_output_file = sample_outdir / "denovo_paths.txt";
    buffer.write_to_file(denovo_output_file);
    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "De novo variant paths written to "
                            << denovo_output_file.string();
}

void pandora_discover_core(const SampleData& sample,
    const std::shared_ptr<Index>& index,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs, const DiscoverOptions& opt)
{
    const auto& sample_name = sample.first;
    const auto& sample_fpath = sample.second;

    // make output dir for this sample
    const auto sample_outdir { opt.outdir / sample_name };
    fs::create_directories(sample_outdir);

    // create kmer graph dir
    const auto kmer_graph_dir { sample_outdir / "kmer_graphs" };
    if (opt.output_kg) {
        fs::create_directories(kmer_graph_dir);
    }

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Constructing pangenome::Graph from read file "
                            << sample_fpath << " (this will take a while)";
    auto pangraph = std::make_shared<pangenome::Graph>();
    uint32_t covg
        = pangraph_from_read_file(sample, pangraph, index, prgs, opt.window_size,
            opt.kmer_size, opt.max_diff, opt.error_rate, opt.min_cluster_size,
            opt.genome_size, opt.illumina, opt.clean, opt.max_covg, opt.threads, sample_outdir);

    const auto pangraph_gfa { sample_outdir / "pandora.pangraph.gfa" };
    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Writing pangenome::Graph to file " << pangraph_gfa;
    write_pangraph_gfa(pangraph_gfa, pangraph);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(warning)
            << "[Sample " << sample_name << "] "
            << "Found no LocalPRGs in the reads for sample " << sample_name;
    }

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Updating local PRGs with hits...";
    pangraph->add_hits_to_kmergraphs();

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Find PRG paths and write to files...";

    // paralell region!
    // shared variable - synced with critical(consensus_fq)
    Fastaq consensus_fq(true, true);

    // shared variable - will denote which nodes we have to remove after the
    // parallel loop. synced with critical(nodes_to_remove)
    std::vector<pangenome::NodePtr> nodes_to_remove;
    nodes_to_remove.reserve(pangraph->nodes.size());

    // transforms the pangraph->nodes from map to vector so that we can run it in
    // parallel
    // TODO: use OMP task instead?
    std::vector<pangenome::NodePtr> pangraphNodesAsVector;
    pangraphNodesAsVector.reserve(pangraph->nodes.size());
    for (auto pan_id_to_node_mapping = pangraph->nodes.begin();
         pan_id_to_node_mapping != pangraph->nodes.end(); ++pan_id_to_node_mapping) {
        pangraphNodesAsVector.push_back(pan_id_to_node_mapping->second);
    }

    std::map<std::string, std::string> locus_to_reads = get_locus_to_reads(
        sample_outdir, sample_name);

    const fs::path denovo_outdir { sample_outdir / "denovo" };
    fs::create_directories(denovo_outdir);

    std::vector<std::string> all_denovo_sequences;
    all_denovo_sequences.reserve(pangraph->nodes.size()*2 + 500);
#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 10)
    for (uint32_t i = 0; i < pangraphNodesAsVector.size(); ++i) {
        // add some progress
        if (i && i % 100 == 0) {
            BOOST_LOG_TRIVIAL(info)
                << "[Sample " << sample_name << "] "
                << ((double)i) / pangraphNodesAsVector.size() * 100 << "% done";
        }

        // get the node
        const auto& pangraph_node = pangraphNodesAsVector[i];

        // add consensus path to fastaq
        std::vector<KmerNodePtr> kmp;
        std::vector<LocalNodePtr> lmp;
        prgs[pangraph_node->prg_id]->add_consensus_path_to_fastaq(consensus_fq,
            pangraph_node, kmp, lmp, opt.window_size, opt.binomial, covg,
            opt.max_num_kmers_to_avg, 0);

        if (kmp.empty()) {
            // mark the node as to remove
#pragma omp critical(nodes_to_remove)
            {
                nodes_to_remove.push_back(pangraph_node);
            }
            continue;
        }

        const std::string& locus = pangraph_node->name;
        const bool no_reads_mapped = locus_to_reads[locus].empty();
        if (no_reads_mapped) {
#pragma omp critical(nodes_to_remove)
            {
                nodes_to_remove.push_back(pangraph_node);
            }
            continue;
        }

        if (opt.output_kg) {
            pangraph_node->kmer_prg_with_coverage.save(
                kmer_graph_dir / (pangraph_node->get_name() + ".kg.gfa"),
                prgs[pangraph_node->prg_id]);
        }

        BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << ", Locus " << locus << "] "
                                << "Running racon to get denovo path";
        const string lmp_seq = prgs[pangraph_node->prg_id]->string_along_path(lmp);

        // builds a mem_fd with the locus reads
        // std::pair<int, std::string> read_locus_fd_and_filepath = build_memfd(locus_to_reads[locus]);
        // TODO: use mem_fd back
        std::string locus_reads_filepath;
        {
            std::stringstream ss;
            ss << locus << ".reads.fa";
            locus_reads_filepath = ss.str();
        }
        locus_reads_filepath = (denovo_outdir / locus_reads_filepath).string();
        build_file(locus_reads_filepath, locus_to_reads[locus]);

        Racon racon(locus, lmp_seq, denovo_outdir, locus_reads_filepath);
        const std::string &polished_sequence = racon.get_polished_sequence();

#pragma omp critical(all_denovo_sequences)
        {
            all_denovo_sequences.push_back(">" + locus);
            all_denovo_sequences.push_back(polished_sequence);
        }
    }

    // remove the nodes marked as to be removed
    for (const auto& node_to_remove : nodes_to_remove) {
        pangraph->remove_node(node_to_remove);
    }

    // write denovo file
    fs::path denovo_filepath = sample_outdir / "denovo_sequences.fa";
    ofstream denovo_filehandler;
    open_file_for_writing(denovo_filepath.string(), denovo_filehandler);
    for (const std::string &line : all_denovo_sequences) {
        denovo_filehandler << line << std::endl;
    }
    denovo_filehandler.close();

    fs::remove_all(denovo_outdir);

    consensus_fq.save(sample_outdir / "pandora.consensus.fq.gz");

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "[Sample " << sample_name << "] "
            << "All nodes which were found have been removed during cleaning. Is "
               "your genome_size accurate?"
            << " Genome size is assumed to be " << opt.genome_size
            << " and can be updated with --genome_size";
    }

    if (opt.output_mapped_read_fa) {
        pangraph->save_mapped_read_strings(sample_fpath, sample_outdir);
    }

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Done discovering!";
}

int pandora_discover(DiscoverOptions& opt)
{
    // checks if minimap2 and racon are available
    if (not tool_exists("minimap2 -h")) {
        throw std::logic_error("minimap2 not found");
    }
    if (not tool_exists("racon -h")) {
        throw std::logic_error("racon not found");
    }

    auto log_level = boost::log::trivial::info;
    if (opt.verbosity == 1) {
        log_level = boost::log::trivial::debug;
    } else if (opt.verbosity > 1) {
        log_level = boost::log::trivial::trace;
    }

    // this is done so that everytime we write a log message, we flush the log
    // if we don't do this, child processes will have messages buffered in their log
    // object when forked and will output repeated log messages
    boost::log::add_console_log(std::cout, boost::log::keywords::auto_flush = true);

    boost::log::core::get()->set_filter(boost::log::trivial::severity >= log_level);

    // =========
    // todo: this all seems strange
    if (opt.error_rate < 0.01) {
        opt.illumina = true;
    }
    if (opt.error_rate > 0.05 and opt.illumina) {
        opt.error_rate = 0.001;
    }
    if (opt.illumina and opt.error_rate > 0.1) {
        opt.error_rate = 0.001;
    }
    if (opt.illumina and opt.max_diff > 200) {
        opt.max_diff = 2 * opt.kmer_size + 1;
    }
    // ==========

    if (opt.window_size > opt.kmer_size) {
        throw std::logic_error("W must NOT be greater than K");
    }
    if (opt.window_size <= 0) {
        throw std::logic_error("W must be a positive integer");
    }
    if (opt.kmer_size <= 0) {
        throw std::logic_error("K must be a positive integer");
    }

    BOOST_LOG_TRIVIAL(info) << "Loading Index and LocalPRGs from file...";
    auto index = std::make_shared<Index>();
    index->load(opt.prgfile, opt.window_size, opt.kmer_size);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);
    load_PRG_kmergraphs(prgs, opt.window_size, opt.kmer_size, opt.prgfile);

    BOOST_LOG_TRIVIAL(info) << "Loading read index file...";
    std::vector<SampleData> samples
        = load_read_index(opt.reads_idx_file);

    // for each sample, run pandora discover
    for (const SampleData& sample : samples) {
        pandora_discover_core(sample, index, prgs, opt);
    }

    // concatenate all denovo files
//    std::vector<fs::path> denovo_paths_files;
//    for (uint32_t sample_id = 0; sample_id < samples.size(); sample_id++) {
//        const auto& sample = samples[sample_id];
//        const auto& sample_name = sample.first;
//        fs::path denovo_output_file = opt.outdir / sample_name / "denovo_paths.txt";
//        denovo_paths_files.push_back(denovo_output_file);
//    }
//    fs::path denovo_output_file = opt.outdir / "denovo_paths.txt";
//    concatenate_denovo_files(denovo_output_file, denovo_paths_files);
//    BOOST_LOG_TRIVIAL(info) << "De novo variant paths written to "
//                            << denovo_output_file.string();
    BOOST_LOG_TRIVIAL(info) << "All done!";

    return 0;
}
