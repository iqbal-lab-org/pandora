//
// Created by michael on 26/9/20.
//

#include "denovo_discovery/discover_main.h"

void setup_discover_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<DiscoverOptions>();
    auto* discover_subcmd = app.add_subcommand("discover",
        "Quasi-map reads to an indexed PRG, infer the "
        "sequence of present loci in the sample and discover novel variants.");

    discover_subcmd
        ->add_option("<TARGET>", opt->prgfile, "An indexed PRG file (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    discover_subcmd
        ->add_option(
            "<QUERY>", opt->readsfile, "Fast{a,q} file containing reads to quasi-map")
        ->required()
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
        ->check(CLI::Range(0, MAX_DENOVO_K))
        ->type_name("INT");

    std::string description = "Max. insertion size for novel variants. Warning: "
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

    description = "Padding either side of candidate variant intervals";
    discover_subcmd->add_option("-P,--pad", opt->candidate_padding, description)
        ->capture_default_str()
        ->type_name("INT");

    description = "Merge candidate variant intervals within distance";
    discover_subcmd->add_option("-d,--merge", opt->merge_dist, description)
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

    fs::create_directories(opt.outdir);
    if (opt.output_kg)
        fs::create_directories(opt.outdir + "/kmer_graphs");

    BOOST_LOG_TRIVIAL(info) << "Loading Index and LocalPRGs from file...";
    auto index = std::make_shared<Index>();
    index->load(opt.prgfile, opt.window_size, opt.kmer_size);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);
    load_PRG_kmergraphs(prgs, opt.window_size, opt.kmer_size, opt.prgfile);

    BOOST_LOG_TRIVIAL(info)
        << "Constructing pangenome::Graph from read file (this will take a while)...";
    auto pangraph = std::make_shared<pangenome::Graph>();
    uint32_t covg
        = pangraph_from_read_file(opt.readsfile, pangraph, index, prgs, opt.window_size,
            opt.kmer_size, opt.max_diff, opt.error_rate, opt.min_cluster_size,
            opt.genome_size, opt.illumina, opt.clean, opt.max_covg, opt.threads);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(info) << "Found none of the LocalPRGs in the reads.";
        BOOST_LOG_TRIVIAL(info) << "Done!";
        return 0;
    }

    BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file " << opt.outdir
                            << "/pandora.pangraph.gfa";
    write_pangraph_gfa(opt.outdir + "/pandora.pangraph.gfa", pangraph);

    BOOST_LOG_TRIVIAL(info) << "Updating local PRGs with hits...";
    pangraph->add_hits_to_kmergraphs(prgs);

    BOOST_LOG_TRIVIAL(info) << "Find PRG paths and write to files...";

    // paralell region!
    // shared variable - synced with critical(consensus_fq)
    Fastaq consensus_fq(true, true);

    // shared variable - synced with critical(candidate_regions)
    CandidateRegions candidate_regions;

    // shared variable - will denote which nodes we have to remove after the
    // parallel loop synced with critical(nodes_to_remove)
    std::vector<pangenome::NodePtr> nodes_to_remove;
    nodes_to_remove.reserve(pangraph->nodes.size());

    Discover discover { opt.min_candidate_covg, opt.min_candidate_len,
        opt.max_candidate_len, opt.candidate_padding, opt.merge_dist };

    // transforms the pangraph->nodes from map to vector so that we can run it in
    // parallel
    // TODO: use OMP task instead?
    std::vector<pangenome::NodePtr> pangraphNodesAsVector;
    pangraphNodesAsVector.reserve(pangraph->nodes.size());
    for (auto pan_id_to_node_mapping = pangraph->nodes.begin();
         pan_id_to_node_mapping != pangraph->nodes.end(); ++pan_id_to_node_mapping) {
        pangraphNodesAsVector.push_back(pan_id_to_node_mapping->second);
    }

    // TODO: check the batch size
#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 10)
    for (uint32_t i = 0; i < pangraphNodesAsVector.size(); ++i) {
        // add some progress
        if (i && i % 100 == 0) {
            BOOST_LOG_TRIVIAL(info)
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
            // pan_id_to_node_mapping = pangraph->remove_node(pangraph_node);
            // mark the node as to remove
#pragma omp critical(nodes_to_remove)
            {
                nodes_to_remove.push_back(pangraph_node);
            }
            continue;
        }

        if (opt.output_kg) {
            pangraph_node->kmer_prg_with_coverage.save(
                opt.outdir + "/kmer_graphs/" + pangraph_node->get_name() + ".kg.gfa",
                prgs[pangraph_node->prg_id]);
        }

        BOOST_LOG_TRIVIAL(info) << "Searching for regions with evidence of novel "
                                   "variants...";
        const TmpPanNode pangraph_node_components { pangraph_node,
            prgs[pangraph_node->prg_id], kmp, lmp };
        auto candidate_regions_for_pan_node {
            discover.find_candidate_regions_for_pan_node(pangraph_node_components)
        };

#pragma omp critical(candidate_regions)
        {
            candidate_regions.insert(candidate_regions_for_pan_node.begin(),
                candidate_regions_for_pan_node.end());
        }
    }

    // build the pileup for candidate regions multithreadly
    BOOST_LOG_TRIVIAL(info) << "Building read pileups for " << candidate_regions.size()
                            << " candidate de novo regions...";
    // the pileup_construction_map function is intentionally left
    // single threaded since it would require too much synchronization
    const auto pileup_construction_map
        = discover.pileup_construction_map(candidate_regions);

    discover.load_candidate_region_pileups(
        opt.readsfile, candidate_regions, pileup_construction_map, opt.threads);

    // remove the nodes marked as to be removed
    for (const auto& node_to_remove : nodes_to_remove) {
        pangraph->remove_node(node_to_remove);
    }

    consensus_fq.save(opt.outdir + "/pandora.consensus.fq.gz");

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "All nodes which were found have been removed during cleaning. Is "
               "your genome_size accurate?"
            << " Genome size is assumed to be " << opt.genome_size
            << " and can be updated with --genome_size";
        return 0;
    }

    DenovoDiscovery denovo { opt.denovo_kmer_size, opt.error_rate,
        opt.max_insertion_size, opt.min_covg_for_node_in_assembly_graph,
        opt.clean_dbg };
    const fs::path denovo_output_directory { fs::path(opt.outdir) / "denovo_paths" };
    fs::create_directories(denovo_output_directory);
    BOOST_LOG_TRIVIAL(info)
        << "Generating de novo variants as paths through their local graph...";

    // TODO: this is hard to parallelize due to GATB's temp files
    for (auto& element : candidate_regions) {
        auto& candidate_region { element.second };
        denovo.find_paths_through_candidate_region(
            candidate_region, denovo_output_directory);
        candidate_region.write_denovo_paths_to_file(denovo_output_directory);
    }
    BOOST_LOG_TRIVIAL(info) << "De novo variant paths written to "
                            << denovo_output_directory.string();

    if (opt.output_mapped_read_fa) {
        pangraph->save_mapped_read_strings(opt.readsfile, opt.outdir);
    }

    BOOST_LOG_TRIVIAL(info) << "Done!";
    return 0;
}
