//
// Created by michael on 26/9/20.
//

#include "denovo_discovery/discover_main.h"
#include "denovo_discovery/racon.h"
#include "denovo_discovery/denovo_record.h"

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
        = pangraph_from_read_file(sample, pangraph, index, prgs,
        opt.window_size, opt.kmer_size, opt.max_diff, opt.error_rate, sample_outdir,
        opt.min_cluster_size, opt.genome_size, opt.illumina, opt.clean, opt.max_covg,
        opt.threads);

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

    std::vector<std::string> all_denovo_sequences;
    all_denovo_sequences.reserve(pangraph->nodes.size()*2 + 500);

    // find denovo paths
    const fs::path denovo_outdir { sample_outdir / "denovo" };
    fs::create_directories(denovo_outdir);

    std::ofstream denovo_paths_out_core_file;
    open_file_for_writing((denovo_outdir / "denovo_paths_core.txt").string(),
        denovo_paths_out_core_file);

    uint32_t number_of_loci_with_denovo_variants = 0;
#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 1) default(shared)
    for (uint32_t i = 0; i < pangraphNodesAsVector.size(); ++i) {
        // add some progress
        if (i && i % 100 == 0) {
            BOOST_LOG_TRIVIAL(info)
                << "[Sample " << sample_name << "] "
                << ((double)i) / pangraphNodesAsVector.size() * 100 << "% done";
        }

        // get the node
        const auto& pangraph_node = pangraphNodesAsVector[i];
        const std::string& locus = pangraph_node->name;

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

        // [TODO RACON]: switch to memfd to improve performance whenever we can
        // builds a mem_fd with the locus reads
        // std::pair<int, std::string> read_locus_fd_and_filepath = build_memfd(locus_to_reads[locus]);

        // build the reads file on disk
        std::string locus_reads_filepath;
        {
            std::stringstream ss;
            ss << locus << ".reads.fa";
            locus_reads_filepath = ss.str();
        }
        locus_reads_filepath = (denovo_outdir / locus_reads_filepath).string();
        build_file(locus_reads_filepath, locus_to_reads[locus]);

        const std::string lmp_seq = prgs[pangraph_node->prg_id]->string_along_path(lmp);
        Racon racon(opt.illumina, opt.kmer_size, locus, lmp_seq,
            denovo_outdir, locus_reads_filepath);
        const std::string &polished_sequence = racon.get_polished_sequence();

        // fs::remove(locus_reads_filepath);

#pragma omp critical(all_denovo_sequences)
        {
            all_denovo_sequences.push_back(">" + locus + " sample=" + sample_name +
                " denovo_sequence");
            all_denovo_sequences.push_back(polished_sequence);
        }

        const std::vector<DenovoVariantRecord> denovo_variants =
            DenovoVariantRecord::get_variants_from_pair_of_sequences(lmp_seq, polished_sequence);

        const bool we_have_denovo_variants = denovo_variants.size() > 0;
        if (we_have_denovo_variants) {
            const std::string lmp_as_string = LocalNode::to_string_vector(lmp);
#pragma omp critical(denovo_paths_out_core_file)
            {
                ++number_of_loci_with_denovo_variants;
                denovo_paths_out_core_file << locus << std::endl;
                denovo_paths_out_core_file << lmp_as_string << std::endl;
                denovo_paths_out_core_file << denovo_variants.size()
                                           << " denovo variants for this locus" << std::endl;
                for (const DenovoVariantRecord &denovo_variant : denovo_variants) {
                    denovo_paths_out_core_file << denovo_variant.to_string() << std::endl;
                }
            }
        }
    }
    denovo_paths_out_core_file.close();

    write_denovo_header_file(sample_name, denovo_outdir,
        number_of_loci_with_denovo_variants);

    fs::path denovo_paths_out_file(sample_outdir / "denovo_paths.txt");
    std::vector<fs::path> denovo_paths_intermediary_files{
        denovo_outdir / "denovo_paths_header.txt",
        denovo_outdir / "denovo_paths_core.txt"
    };
    concatenate_text_files(denovo_paths_out_file, denovo_paths_intermediary_files);
    
    // remove the nodes marked as to be removed
    for (const auto& node_to_remove : nodes_to_remove) {
        pangraph->remove_node(node_to_remove);
    }

    // write denovo file
    fs::path denovo_filepath = sample_outdir / "denovo_sequences.fa";
    std::ofstream denovo_filehandler;
    open_file_for_writing(denovo_filepath.string(), denovo_filehandler);
    for (const std::string &line : all_denovo_sequences) {
        denovo_filehandler << line << std::endl;
    }
    denovo_filehandler.close();

    // fs::remove_all(denovo_outdir);

    consensus_fq.save(sample_outdir / "pandora.consensus.fq.gz");

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "[Sample " << sample_name << "] "
            << "All nodes which were found have been removed during cleaning. Is "
               "your genome_size accurate?"
            << " Genome size is assumed to be " << opt.genome_size
            << " and can be updated with --genome_size";
    }

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Done discovering!";
}

int pandora_discover(DiscoverOptions& opt)
{
    auto log_level = boost::log::trivial::info;
    if (opt.verbosity == 1) {
        log_level = boost::log::trivial::debug;
    } else if (opt.verbosity > 1) {
        log_level = boost::log::trivial::trace;
    }
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

    concatenate_all_denovo_files(samples, opt.outdir);
    BOOST_LOG_TRIVIAL(info) << "All done!";

    return 0;
}
