/*
 * Todo:
 * 1. Allow reverse complement strand minimizers
 * 2. Make sure mapped region not got large indels from unsupported intervals
 * 3. Command line argument handling for option of dumped_db screening
 * 4. Change structure so not copying minimzers etc, but use pointers
 */
#include "map_main.h"

void setup_map_subcommand(CLI::App& app)
{
    // todo: make groups for opts
    auto opt = std::make_shared<MapOptions>();
    auto map_subcmd = app.add_subcommand("map",
        "Quasi-map reads to an indexed PRG, infer the "
        "sequence of present loci in the sample, and (optionally) genotype/discover "
        "variants.");

    map_subcmd
        ->add_option("<TARGET>", opt->prgfile, "An indexed PRG file (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    map_subcmd
        ->add_option(
            "<QUERY>", opt->readsfile, "Fast{a,q} file containing reads to quasi-map")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    std::stringstream desc;
    desc << "Window size for (w,k)-minimizers (must be <=k) [default: "
         << opt->window_size << "]";
    map_subcmd->add_option("-w", opt->window_size, desc.str())
        ->type_name("INT")
        ->group("Indexing");

    desc.str(std::string());
    desc << "K-mer size for (w,k)-minimizers [default: " << opt->kmer_size << "]";
    map_subcmd->add_option("-k", opt->kmer_size, desc.str())
        ->type_name("INT")
        ->group("Indexing");

    desc.str(std::string());
    desc << "Directory to write output files to [default: " << opt->outdir << "]";
    map_subcmd->add_option("-o,--outdir", opt->outdir, desc.str())
        ->type_name("DIR")
        ->group("Input/Output");

    desc.str(std::string());
    desc << "Maximum number of threads to use [default: " << opt->threads << "]";
    map_subcmd->add_option("-t,--threads", opt->threads, desc.str())
        ->type_name("INT")
        ->group("Input/Output");

    desc.str(std::string());
    desc << "Fasta file with an entry for each loci. The entries will be used as the "
            "reference sequence for their respective loci. The sequence MUST have a "
            "perfect match in <TARGET> and the same name";
    map_subcmd->add_option("--vcf-refs", opt->vcf_refs_file, desc.str())
        ->type_name("FILE")
        ->group("Input/Output");

    desc.str(std::string());
    desc << "Estimated error rate for reads [default: " << opt->error_rate << "]";
    map_subcmd->add_option("-e,--error-rate", opt->error_rate, desc.str())
        ->group("Parameter Estimation");

    // todo: how necessary is this opt if we remove max_cog?
    desc.str(std::string());
    desc << "Estimated length of the genome - used for coverage estimation [default: "
         << opt->genome_size << "]";
    map_subcmd->add_option("-g,--genome-size", opt->genome_size, desc.str())
        ->type_name("INT")
        ->group("Parameter Estimation");

    desc.str(std::string());
    desc << "Maximum distance (bp) between consecutive hits within a cluster [default: "
         << opt->max_diff << "]";
    map_subcmd->add_option("-m,--max-diff", opt->max_diff, desc.str())
        ->type_name("INT")
        ->group("Mapping");

    // todo: suggest rename
    map_subcmd
        ->add_flag("--output-kg", opt->output_kg,
            "Save kmer graphs with forward and reverse coverage annotations for found "
            "loci")
        ->group("Input/Output");

    // todo: suggest rename
    map_subcmd
        ->add_flag(
            "--output-vcf", opt->output_vcf, "Save a VCF file for each found loci")
        ->group("Input/Output");

    // todo: suggest rename
    map_subcmd
        ->add_flag("--output-comparison-paths", opt->output_comparison_paths,
            "Save a fasta file for a random selection of paths through loci")
        ->group("Input/Output");

    // todo: suggest rename
    map_subcmd
        ->add_flag("--output-covgs", opt->output_covgs,
            "Save a file of coverages for each loci present - one number per base")
        ->group("Input/Output");

    // todo: suggest rename
    map_subcmd
        ->add_flag("--output-mapped-read-fa", opt->output_mapped_read_fa,
            "Save a file for each loci containing read parts which overlapped it")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("--illumina", opt->illumina,
            "Reads are from Illumina. Alters error rate used and adjusts for shorter "
            "reads")
        ->group("Preset");

    map_subcmd
        ->add_flag(
            "--clean", opt->clean, "Add a step to clean and detangle the pangraph")
        ->group("Filtering");

    map_subcmd
        ->add_flag("--bin", opt->binomial,
            "Use binomial model for kmer coverages [default: negative binomial]")
        ->group("Parameter Estimation");

    // todo: suggest removing this parameter
    desc.str(std::string());
    desc << "Maximum average coverage of reads to accept [default: " << opt->max_covg
         << "]";
    map_subcmd->add_option("--max-covg", opt->max_covg, desc.str())
        ->type_name("INT")
        ->group("Filtering");

    auto genotype_validator = [](const std::string& str) {
        bool valid_genotype_option = str == "global" || str == "local";
        if (!valid_genotype_option) {
            return std::string("--genotype should be either 'global' or 'local'.");
        }
        return std::string();
    };
    // todo: zam mentioned we probably only want global - remove local and make this a
    // flag?
    map_subcmd
        ->add_option("--genotype", opt->genotype,
            "Add extra step to carefully genotype sites. There are two modes: 'global' "
            "(ML path oriented) or 'local' (coverage oriented)")
        ->type_name("MODE")
        ->check(genotype_validator)
        ->group("Consensus/Variant Calling");

    map_subcmd
        ->add_flag("--snps", opt->snps_only, "When genotyping, only include SNP sites")
        ->group("Consensus/Variant Calling");

    map_subcmd
        ->add_flag(
            "-d,--discover", opt->discover, "Add a step to discover de novo variants")
        ->group("Consensus/Variant Calling");

    desc.str(std::string());
    desc << "K-mer size to use when disovering de novo variants [default: "
         << std::to_string(opt->denovo_kmer_size) << "]";
    map_subcmd->add_option("--discover-k", opt->denovo_kmer_size, desc.str())
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    desc.str(std::string());
    desc << "Minimum size of a cluster of hits between a read and a loci to consider "
            "the loci present [default: "
         << opt->min_cluster_size << "]";
    map_subcmd->add_option("-c,--min-cluster-size", opt->min_cluster_size, desc.str())
        ->type_name("INT")
        ->group("Mapping");

    desc.str(std::string());
    desc << "Maximum number of kmers to average over when selecting the maximum "
            "likelihood path [default: "
         << opt->max_num_kmers_to_avg << "]";
    map_subcmd->add_option("--kmer-avg", opt->max_num_kmers_to_avg, desc.str())
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    desc.str(std::string());
    desc << "Should this be exposed?"; // todo
    map_subcmd->add_option("--min-kmer-covg", opt->min_kmer_covg, desc.str())
        ->type_name("INT")
        ->group("??");

    desc.str(std::string());
    desc << "Should this be exposed?"; // todo
    map_subcmd->add_option("--min-allele-covg-gt", opt->min_allele_covg_gt, desc.str())
        ->type_name("INT")
        ->group("??");

    desc.str(std::string());
    desc << "Should this be exposed?"; // todo
    map_subcmd->add_option("--min-total-covg-gt", opt->min_total_covg_gt, desc.str())
        ->type_name("INT")
        ->group("??");

    desc.str(std::string());
    desc << "Should this be exposed?"; // todo
    map_subcmd->add_option("--min-diff-covg-gt", opt->min_diff_covg_gt, desc.str())
        ->type_name("INT")
        ->group("??");

    desc.str(std::string());
    desc << "Should this be exposed?"; // todo
    map_subcmd
        ->add_option("--min-allele-fraction-covg-gt", opt->min_allele_fraction_covg_gt,
            desc.str())
        ->type_name("INT")
        ->group("??");

    desc.str(std::string());
    desc << ""; // todo
    map_subcmd
        ->add_option("--genotyping-error-rate", opt->genotyping_error_rate, desc.str())
        ->group("Consensus/Variant Calling");

    desc.str(std::string());
    desc << ""; // todo
    map_subcmd
        ->add_option("--confidence-threshold", opt->confidence_threshold, desc.str())
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    map_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    map_subcmd->callback([opt]() { pandora_map(*opt); });
}

int pandora_map(MapOptions& opt)
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

    bool do_global_genotyping { opt.genotype == "global" };
    bool do_local_genotyping { opt.genotype == "local" };

    if (do_global_genotyping or do_local_genotyping)
        opt.output_vcf = true;

    GenotypingOptions genotyping_options({}, opt.genotyping_error_rate,
        opt.confidence_threshold, opt.min_allele_covg_gt,
        opt.min_allele_fraction_covg_gt, opt.min_total_covg_gt, opt.min_diff_covg_gt, 0,
        false);

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
        BOOST_LOG_TRIVIAL(info) << "Found non of the LocalPRGs in the reads.";
        BOOST_LOG_TRIVIAL(info) << "Done!";
        return 0;
    }

    BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file " << opt.outdir
                            << "/pandora.pangraph.gfa";
    write_pangraph_gfa(opt.outdir + "/pandora.pangraph.gfa", pangraph);

    BOOST_LOG_TRIVIAL(info) << "Update LocalPRGs with hits...";
    uint32_t sample_id = 0;
    pangraph->add_hits_to_kmergraphs(prgs);

    BOOST_LOG_TRIVIAL(info) << "Estimate parameters for kmer graph model...";
    auto exp_depth_covg = estimate_parameters(pangraph, opt.outdir, opt.kmer_size,
        opt.error_rate, covg, opt.binomial, sample_id);
    genotyping_options.add_exp_depth_covg(exp_depth_covg);
    if (genotyping_options.get_min_kmer_covg() == 0) {
        genotyping_options.set_min_kmer_covg(exp_depth_covg / 10);
    }

    BOOST_LOG_TRIVIAL(info) << "Find PRG paths and write to files...";

    // paralell region!
    // shared variable - synced with critical(consensus_fq)
    Fastaq consensus_fq(true, true);

    // shared variable - synced with critical(master_vcf)
    VCF master_vcf(&genotyping_options);

    // shared variable - synced with critical(candidate_regions)
    CandidateRegions candidate_regions;

    // shared variable - will denote which nodes we have to remove after the
    // parallel loop synced with critical(nodes_to_remove)
    std::vector<pangenome::NodePtr> nodes_to_remove;
    nodes_to_remove.reserve(pangraph->nodes.size());

    // this a read-only var, no need for sync
    VCFRefs vcf_refs;
    if (opt.output_vcf and !opt.vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(opt.vcf_refs_file, vcf_refs);
    }

    // transforms the pangraph->nodes from map to vector so that we can run it in
    // parallel
    // TODO: use OMP task instead?
    std::vector<pangenome::NodePtr> pangraphNodesAsVector;
    pangraphNodesAsVector.reserve(pangraph->nodes.size());
    for (auto pan_id_to_node_mapping = pangraph->nodes.begin();
         pan_id_to_node_mapping != pangraph->nodes.end(); ++pan_id_to_node_mapping)
        pangraphNodesAsVector.push_back(pan_id_to_node_mapping->second);

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

        // get the vcf_ref, if applicable
        std::string vcf_ref;
        if (opt.output_vcf and !opt.vcf_refs_file.empty()
            and vcf_refs.find(prgs[pangraph_node->prg_id]->name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prgs[pangraph_node->prg_id]->name];
        }

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

        if (opt.output_vcf) {
            // TODO: this takes a lot of time and should be optimized, but it is
            // only called in this part, so maybe this should be low prioritized
            prgs[pangraph_node->prg_id]->add_variants_to_vcf(
                master_vcf, pangraph_node, vcf_ref, kmp, lmp, opt.min_kmer_covg);
        }

        if (opt.discover) {
            BOOST_LOG_TRIVIAL(info) << "Searching for regions with evidence of de novo "
                                       "variants...";
            const TmpPanNode pangraph_node_components { pangraph_node,
                prgs[pangraph_node->prg_id], kmp, lmp };
            auto candidate_regions_for_pan_node { find_candidate_regions_for_pan_node(
                pangraph_node_components, opt.denovo_kmer_size * 2) };

#pragma omp critical(candidate_regions)
            {
                candidate_regions.insert(candidate_regions_for_pan_node.begin(),
                    candidate_regions_for_pan_node.end());
            }
        }
    }

    // build the pileup for candidate regions multithreadly
    if (opt.discover) {
        BOOST_LOG_TRIVIAL(info)
            << "Building read pileups for candidate de novo regions...";
        // the construct_pileup_construction_map function is intentionally left
        // single threaded since it would require too much synchronization
        const auto pileup_construction_map
            = construct_pileup_construction_map(candidate_regions);

        load_all_candidate_regions_pileups_from_fastq(
            opt.readsfile, candidate_regions, pileup_construction_map, opt.threads);
    }

    // remove the nodes marked as to be removed
    for (const auto& node_to_remove : nodes_to_remove)
        pangraph->remove_node(node_to_remove);

    consensus_fq.save(opt.outdir + "/pandora.consensus.fq.gz");
    if (opt.output_vcf)
        master_vcf.save(opt.outdir + "/pandora_consensus.vcf", true, false);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "All nodes which were found have been removed during cleaning. Is "
               "your genome_size accurate?"
            << " Genome size is assumed to be " << opt.genome_size
            << " and can be updated with --genome_size";
        return 0;
    }

    if (do_global_genotyping or do_local_genotyping) {
        master_vcf.genotype(do_global_genotyping, do_local_genotyping);
        if (opt.snps_only)
            master_vcf.save(opt.outdir + "/pandora_genotyped_" + opt.genotype + ".vcf",
                false, true, false, true, true, true, true, false, false, false);
        else
            master_vcf.save(opt.outdir + "/pandora_genotyped_" + opt.genotype + ".vcf",
                false, true);
    }

    if (opt.discover) {
        DenovoDiscovery denovo { opt.denovo_kmer_size, opt.error_rate };
        const fs::path denovo_output_directory { fs::path(opt.outdir)
            / "denovo_paths" };
        fs::create_directories(denovo_output_directory);
        BOOST_LOG_TRIVIAL(info)
            << "Generating de novo variants as paths through their local graph...";

        for (auto& element : candidate_regions) {
            auto& candidate_region { element.second };
            denovo.find_paths_through_candidate_region(candidate_region,
                denovo_output_directory); // TODO: this is hard to parallelize due
                                          // to GATB's temp files
            candidate_region.write_denovo_paths_to_file(denovo_output_directory);
        }
        BOOST_LOG_TRIVIAL(info)
            << "De novo variant paths written to " << denovo_output_directory.string();
    }

    if (opt.output_mapped_read_fa) {
        pangraph->save_mapped_read_strings(opt.readsfile, opt.outdir);
    }

    BOOST_LOG_TRIVIAL(info) << "Done!";
    return 0;
}
