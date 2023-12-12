#include "map_main.h"
#include "cli_helpers.h"

void setup_map_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<MapOptions>();
    auto* map_subcmd = app.add_subcommand("map",
        "Quasi-map reads to a pandora index, infer the sequence of present loci in the "
        "sample, and optionally genotype variants.");

    map_subcmd
        ->add_option("<TARGET>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
        ->type_name("FILE");

    map_subcmd
        ->add_option(
            "<QUERY>", opt->readsfile, "Fast{a,q} file containing reads to quasi-map")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    map_subcmd
        ->add_option("-o,--outdir", opt->outdir, "Directory to write output files to")
        ->type_name("DIR")
        ->capture_default_str()
        ->transform(make_absolute)
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
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->group("Input/Output");

    map_subcmd
        ->add_option(
            "-e,--error-rate", opt->error_rate, "Estimated error rate for reads")
        ->capture_default_str()
        ->group("Parameter Estimation");

    map_subcmd
        ->add_option("-g,--genome-size", opt->genome_size,
            "Estimated length of the genome - used for coverage estimation. Can pass "
            "string such as 4.4m, 100k etc.")
        ->transform(transform_cli_gsize)
        ->capture_default_str()
        ->type_name("STR/INT")
        ->group("Parameter Estimation");

    map_subcmd
        ->add_option("-m,--max-diff", opt->max_diff,
            "Maximum distance (bp) between consecutive hits within a cluster")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    description
        = "When two clusters of hits are conflicting, the one with highest number of unique minimisers "
          "will be kept. However, if the difference between the number of unique minimisers is too small, "
          "less than this parameter, then we will prefer the cluster that has higher target coverage.";
    map_subcmd
        ->add_option("--cluster-mini-tolerance", opt->conflicting_clusters_minimiser_tolerance, description)
        ->capture_default_str()
        ->type_name("FLOAT")
        ->group("Mapping");

    description
        = "Minimum proportion of overlap between two clusters of hits to consider "
          "them conflicting. Only one of the conflicting clusters will be kept.";
    map_subcmd
        ->add_option("--min-cluster-overlap", opt->conflicting_clusters_overlap_threshold, description)
        ->capture_default_str()
        ->type_name("FLOAT")
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
        ->add_flag("-I,--illumina", opt->illumina,
            "Reads are from Illumina. Alters error rate used and adjusts for shorter "
            "reads")
        ->group("Preset");

    map_subcmd
        ->add_flag("--bin", opt->binomial,
            "Use binomial model for kmer coverages [default: negative binomial]")
        ->group("Parameter Estimation");

    map_subcmd
        ->add_flag("--dont-auto-update-params", opt->do_not_auto_update_params,
            "By default, pandora automatically updates error rate and kmer coverage model parameters based "
            "on the mapping of the previous sample. This could potentially generate "
            "more accurate results if your samples have no sequencing issues and a "
            "consistent protocol was followed for the sequencing of all samples. If this is "
            "not the case, deactivate this feature by activating this flag")
        ->group("Parameter Estimation");

    map_subcmd
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Filtering");

    map_subcmd
        ->add_option(
            "--min-abs-gene-coverage", opt->min_absolute_gene_coverage,
            "Minimum absolute mean gene coverage to keep a gene. Given the "
            "coverage on the kmers of the maximum likelihood path of a gene, we compute "
            "the mean gene coverage and compare with the value in this "
            "parameter. If the mean is lower than this parameter, "
            "the gene is filtered out, e.g. if this parameter value is "
            "3, then all genes with mean <3 will be filtered out.")
        ->capture_default_str()
        ->type_name("FLOAT")
        ->group("Filtering");

    map_subcmd
        ->add_option(
            "--min-rel-gene-coverage", opt->min_relative_gene_coverage,
            "Minimum relative mean gene coverage to keep a gene. This is a proportion, between 0.0 and 1.0. "
            "Given the coverage on the kmers of the maximum likelihood path of a gene, we compute "
            "the mean gene coverage and compare with the value in this "
            "parameter and the global coverage. If the mean is lower"
            " than the computed value, the gene is filtered out, e.g. if this parameter value is "
            "0.05, then all genes with mean < 5% of the global coverage will be "
            "filtered out.")
        ->capture_default_str()
        ->type_name("FLOAT")
        ->group("Filtering");

    map_subcmd
        ->add_option(
            "--max-rel-gene-coverage", opt->max_relative_gene_coverage,
            "Maximum relative mean gene coverage to keep a gene. "
            "Given the coverage on the kmers of the maximum likelihood path of a gene, we compute "
            "the mean gene coverage and compare with the value in this "
            "parameter and the global coverage. If the mean is higher"
            " than the computed value, the gene is filtered out, e.g. if this parameter value is "
            "10, then all genes with mean >10 times the global coverage will be "
            "filtered out.")
        ->capture_default_str()
        ->type_name("FLOAT")
        ->group("Filtering");

    map_subcmd
        ->add_flag(
            "--no-gene-coverage-filtering", opt->no_gene_coverage_filtering,
            "Do not filter genes based on their coverage, effectively ignoring the three "
            "previous params. This is useful if you are not using read datasets.")
        ->group("Filtering");

    map_subcmd
        ->add_option(
            "--min-gene-coverage-proportion", opt->min_gene_coverage_proportion,
            "Minimum gene coverage proportion to keep a gene. "
            "Given the coverage on the kmers of the maximum likelihood path of a gene, we compute "
            "the number of bases that have at least one read covering it. "
            "If the proportion of such bases is larger than the value in this "
            "parameter, the gene is kept. Otherwise, the gene is filtered out.")
        ->capture_default_str()
        ->type_name("FLOAT")
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

    map_subcmd
        ->add_flag("-K,--debugging-files", opt->keep_extra_debugging_files,
            "Keep extra debugging files. Warning: this might "
            "create thousands of files.")
        ->group("Debugging");

    map_subcmd
        ->add_option(
            "-r,--rng-seed", opt->rng_seed, "RNG seed, an int>0 to force deterministic "
                                            "mapping when multiple optimal mappings are "
                                            "possible. To be avoided except in "
                                            "debugging/investigation scenarios. A value "
                                            "of 0 will be interpreted as no seed given "
                                            "and mapping will not be deterministic.")
        ->capture_default_str()
        ->group("Debugging");

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

    if (opt.illumina) {
        opt.error_rate = 0.001;
    }

    if (opt.genotype) {
        opt.output_vcf = true;
    }

    GenotypingOptions genotyping_options({}, opt.genotyping_error_rate,
        opt.confidence_threshold, opt.min_allele_covg_gt,
        opt.min_allele_fraction_covg_gt, opt.min_total_covg_gt, opt.min_diff_covg_gt, 0,
        false);

    BOOST_LOG_TRIVIAL(info) << "Loading Index...";
    Index index = Index::load(opt.index_file.string());
    BOOST_LOG_TRIVIAL(info) << "Index loaded successfully!";

    if (opt.illumina and opt.max_diff > 200) {
        opt.max_diff = 2 * index.get_kmer_size() + 1;
    }

    fs::create_directories(opt.outdir);
    const auto kmer_graphs_dir { opt.outdir / "kmer_graphs" };
    if (opt.output_kg) {
        fs::create_directories(kmer_graphs_dir);
    }

    BOOST_LOG_TRIVIAL(info)
        << "Constructing pangenome::Graph from read file (this will take a while)...";
    const SampleData sample(opt.readsfile.stem().string(), opt.readsfile.string());
    auto pangraph = std::make_shared<pangenome::Graph>();
    uint32_t covg
        = pangraph_from_read_file(sample, pangraph, index, opt.max_diff, opt.error_rate,
            opt.outdir, opt.min_cluster_size, opt.genome_size, opt.max_covg,
            opt.conflicting_clusters_overlap_threshold, opt.conflicting_clusters_minimiser_tolerance,
            opt.threads, opt.keep_extra_debugging_files, opt.rng_seed);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(info) << "Found none of the LocalPRGs in the reads.";
        BOOST_LOG_TRIVIAL(info) << "Done!";
        return 0;
    }

    BOOST_LOG_TRIVIAL(info) << "Estimating parameters for kmer graph model...";
    auto exp_depth_covg = estimate_parameters(pangraph, opt.outdir, index.get_kmer_size(),
        opt.error_rate, covg, opt.binomial, 0, opt.do_not_auto_update_params);
    genotyping_options.add_exp_depth_covg(exp_depth_covg);

    if (genotyping_options.get_min_kmer_covg() == 0) {
        genotyping_options.set_min_kmer_covg(exp_depth_covg / 10);
    }

    BOOST_LOG_TRIVIAL(info) << "Find PRG paths and write to files...";

    // paralell region!
    // shared variable - synced with critical(consensus_fq)
    Fastaq consensus_fq(true, true);

    // shared variable - will denote which nodes we have to remove after the
    // parallel loop synced with critical(nodes_to_remove)
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

    // shared variable - synced with critical(master_vcf)
    VCF master_vcf(&genotyping_options);

    // this a read-only var, no need for sync
    VCFRefs vcf_refs;
    if (opt.output_vcf and !opt.vcf_refs_file.empty()) {
        vcf_refs.reserve(index.get_number_of_prgs());
        load_vcf_refs_file(opt.vcf_refs_file, vcf_refs);
    }

#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 10) default(shared)
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
        const std::string &prg_name = index.get_prg_name_given_id(pangraph_node->prg_id);
        if (opt.output_vcf and !opt.vcf_refs_file.empty()
            and vcf_refs.find(prg_name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prg_name];
        }

        // add consensus path to fastaq
        const auto &prg = index.get_prg_given_id(pangraph_node->prg_id);

        std::vector<KmerNodePtr> kmp;
        std::vector<LocalNodePtr> lmp;
        prg->add_consensus_path_to_fastaq(consensus_fq, pangraph_node, kmp, lmp,
            index.get_window_size(), opt.binomial, covg, opt.max_num_kmers_to_avg, 0,
            opt.min_absolute_gene_coverage, opt.min_relative_gene_coverage,
            opt.max_relative_gene_coverage, opt.no_gene_coverage_filtering);
            opt.max_relative_gene_coverage, opt.min_gene_coverage_proportion);

        if (kmp.empty()) {
#pragma omp critical(nodes_to_remove)
            {
                nodes_to_remove.push_back(pangraph_node);
            }
            continue;
        }

        if (opt.output_kg) {
            pangraph_node->kmer_prg_with_coverage.save(kmer_graphs_dir / (pangraph_node->get_name() + ".kg.gfa"), prg);
        }

        if (opt.output_vcf) {
            // TODO: this takes a lot of time and should be optimized, but it is
            // only called in this part, so maybe this should be low prioritized
            prg->add_variants_to_vcf(master_vcf, pangraph_node, vcf_ref, kmp, lmp);
        }
    }

    // remove the nodes marked as to be removed
    for (const auto& node_to_remove : nodes_to_remove)
        pangraph->remove_node(node_to_remove);

    consensus_fq.save(opt.outdir / "pandora.consensus.fq.gz");
    if (opt.output_vcf) {
        master_vcf.save(opt.outdir / "pandora_consensus.vcf", true, false);
    }

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

        // save the genotyped vcf to disk
        const auto gt_vcf_filepath { opt.outdir / "pandora_genotyped.vcf" };
        if (opt.snps_only) {
            master_vcf.save(gt_vcf_filepath, false, true, false, true, true, true, true,
                false, false, false);
        } else {
            master_vcf.save(gt_vcf_filepath, false, true);
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Done!";
    return 0;
}
