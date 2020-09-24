#include "map_main.h"

void setup_map_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<MapOptions>();
    auto* map_subcmd = app.add_subcommand("map",
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

    map_subcmd
        ->add_option(
            "-w", opt->window_size, "Window size for (w,k)-minimizers (must be <=k)")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    map_subcmd->add_option("-k", opt->kmer_size, "K-mer size for (w,k)-minimizers")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    map_subcmd
        ->add_option("-o,--outdir", opt->outdir, "Directory to write output files to")
        ->type_name("DIR")
        ->capture_default_str()
        ->group("Input/Output");

    map_subcmd
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Input/Output");

    std::string description = "Fasta file with a reference sequence to use for each "
                              "loci. The sequence MUST have a "
                              "perfect match in <TARGET> and the same name";
    map_subcmd->add_option("--vcf-refs", opt->vcf_refs_file, description)
        ->type_name("FILE")
        ->group("Input/Output");

    map_subcmd
        ->add_option(
            "-e,--error-rate", opt->error_rate, "Estimated error rate for reads")
        ->capture_default_str()
        ->group("Parameter Estimation");

    // todo: how necessary is this opt if we remove max_cog?
    map_subcmd
        ->add_option("-g,--genome-size", opt->genome_size,
            "Estimated length of the genome - used for coverage estimation")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Parameter Estimation");

    map_subcmd
        ->add_option("-m,--max-diff", opt->max_diff,
            "Maximum distance (bp) between consecutive hits within a cluster")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    map_subcmd
        ->add_flag("--kg", opt->output_kg,
            "Save kmer graphs with forward and reverse coverage annotations for found "
            "loci")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("--loci-vcf", opt->output_vcf, "Save a VCF file for each found loci")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("-C,--comparison-paths", opt->output_comparison_paths,
            "Save a fasta file for a random selection of paths through loci")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("--coverages", opt->output_covgs,
            "Save a file of coverages for each loci present - one number per base")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("-M,--mapped-reads", opt->output_mapped_read_fa,
            "Save a fasta file for each loci containing read parts which overlapped it")
        ->group("Input/Output");

    map_subcmd
        ->add_flag("-I,--illumina", opt->illumina,
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

    map_subcmd
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Filtering");

    description = "Add extra step to carefully genotype sites.";
    auto* gt_opt = map_subcmd->add_flag("--genotype", opt->genotype, description)
                       ->group("Consensus/Variant Calling");

    description = "(Intended for developers) Use coverage-oriented (local) genotyping "
                  "instead of the default ML path-oriented (global) approach.";
    map_subcmd->add_flag("--local", opt->local_genotype, description)
        ->needs(gt_opt)
        ->group("Genotyping");

    map_subcmd
        ->add_flag("--snps", opt->snps_only, "When genotyping, only include SNP sites")
        ->group("Consensus/Variant Calling");

    map_subcmd
        ->add_flag(
            "-d,--discover", opt->discover, "Add a step to discover de novo variants")
        ->group("Consensus/Variant Calling");

    map_subcmd
        ->add_option("--discover-k", opt->denovo_kmer_size,
            "Kmer size to use when disovering de novo variants")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    description = "Max. insertion size when discovering de novo variants. Warning: "
                  "setting too long may impair performance";
    map_subcmd->add_option("--max-ins", opt->max_insertion_size, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    description = "Minimum node/kmer depth in the de Bruijn graph used for discovering "
                  "de novo variants";
    map_subcmd
        ->add_option(
            "--min-dbg-dp", opt->min_covg_for_node_in_assembly_graph, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    description
        = "Minimum size of a cluster of hits between a read and a loci to consider "
          "the loci present";
    map_subcmd->add_option("-c,--min-cluster-size", opt->min_cluster_size, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    description = "Maximum number of kmers to average over when selecting the maximum "
                  "likelihood path";
    map_subcmd->add_option("--kmer-avg", opt->max_num_kmers_to_avg, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    description
        = "Hard threshold for the minimum allele coverage allowed when genotyping";
    map_subcmd->add_option("-a", opt->min_allele_covg_gt, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

    description = "The minimum required total coverage for a site when genotyping";
    map_subcmd->add_option("-s", opt->min_total_covg_gt, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

    description = "Minimum difference in coverage on a site required between the first "
                  "and second maximum likelihood path";
    map_subcmd->add_option("-D", opt->min_diff_covg_gt, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Genotyping");

    description = "Minimum allele coverage, as a fraction of the expected coverage, "
                  "allowed when genotyping";
    map_subcmd->add_option("-F", opt->min_allele_fraction_covg_gt, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Genotyping");

    description = "When genotyping, assume that coverage on alternative alleles arises "
                  "as a result of an error process with rate -E.";
    map_subcmd
        ->add_option("-E,--gt-error-rate", opt->genotyping_error_rate, description)
        ->capture_default_str()
        ->group("Genotyping");

    description = "Minimum genotype confidence (GT_CONF) required to make a call";
    map_subcmd->add_option("-G,--gt-conf", opt->confidence_threshold, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

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

    if (opt.genotype) {
        opt.output_vcf = true;
    }

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

    BOOST_LOG_TRIVIAL(info) << "Updating local PRGs with hits...";
    uint32_t sample_id = 0;
    pangraph->add_hits_to_kmergraphs(prgs);

    BOOST_LOG_TRIVIAL(info) << "Estimating parameters for kmer graph model...";
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
                master_vcf, pangraph_node, vcf_ref, kmp, lmp);
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
            << "Building read pileups for " << candidate_regions.size()
            << " candidate de novo regions...";
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

    if (opt.genotype) {
        BOOST_LOG_TRIVIAL(info) << "Genotyping VCF...";
        master_vcf.genotype(opt.local_genotype);
        BOOST_LOG_TRIVIAL(info) << "Finished genotyping VCF";
        if (opt.snps_only) {
            master_vcf.save(opt.outdir + "/pandora_genotyped.vcf", false, true, false,
                true, true, true, true, false, false, false);
        } else {
            master_vcf.save(opt.outdir + "/pandora_genotyped.vcf", false, true);
        }
    }

    if (opt.discover) {
        DenovoDiscovery denovo { opt.denovo_kmer_size, opt.error_rate,
            opt.max_insertion_size, opt.min_covg_for_node_in_assembly_graph };
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
