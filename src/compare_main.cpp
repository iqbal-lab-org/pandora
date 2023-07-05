#include "compare_main.h"
#include "cli_helpers.h"

void setup_compare_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<CompareOptions>();
    auto* compare_subcmd = app.add_subcommand("compare",
        "Quasi-map reads from multiple samples to a pandora index, infer the "
        "sequence of present loci in each sample, and call variants between the "
        "samples.");

    compare_subcmd
        ->add_option("<TARGET>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
        ->type_name("FILE");

    std::string description
        = "A tab-delimited file where each line is a sample identifier followed by "
          "the path to the fast{a,q} of reads for that sample";
    compare_subcmd->add_option("<QUERY_IDX>", opt->reads_idx_file, description)
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    compare_subcmd
        ->add_option("-o,--outdir", opt->outdir, "Directory to write output files to")
        ->type_name("DIR")
        ->capture_default_str()
        ->transform(make_absolute)
        ->group("Input/Output");

    compare_subcmd
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Input/Output");

    description = "Fasta file with a reference sequence to use for each loci. The "
                  "sequence MUST have a perfect match in <TARGET> and the same name";
    compare_subcmd->add_option("--vcf-refs", opt->vcf_refs_file, description)
        ->type_name("FILE")
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->group("Input/Output");

    compare_subcmd
        ->add_option(
            "-e,--error-rate", opt->error_rate, "Estimated error rate for reads")
        ->capture_default_str()
        ->group("Parameter Estimation");

    compare_subcmd
        ->add_option("-g,--genome-size", opt->genome_size,
            "Estimated length of the genome - used for coverage estimation. Can pass "
            "string such as 4.4m, 100k etc.")
        ->transform(transform_cli_gsize)
        ->capture_default_str()
        ->type_name("STR/INT")
        ->group("Parameter Estimation");

    compare_subcmd
        ->add_option("-m,--max-diff", opt->max_diff,
            "Maximum distance (bp) between consecutive hits within a cluster")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    compare_subcmd
        ->add_flag("--loci-vcf", opt->output_vcf, "Save a VCF file for each found loci")
        ->group("Input/Output");

    compare_subcmd
        ->add_flag("-I,--illumina", opt->illumina,
            "Reads are from Illumina. Alters error rate used and adjusts for shorter "
            "reads")
        ->group("Preset");

    compare_subcmd
        ->add_flag(
            "--clean", opt->clean, "Add a step to clean and detangle the pangraph")
        ->group("Filtering");

    compare_subcmd
        ->add_flag("--bin", opt->binomial,
            "Use binomial model for kmer coverages [default: negative binomial]")
        ->group("Parameter Estimation");

    compare_subcmd
        ->add_flag("--dont-auto-update-params", opt->do_not_auto_update_params,
            "By default, pandora automatically updates error rate and kmer coverage model parameters based "
            "on the mapping of the previous sample. This could potentially generate "
            "more accurate results if your samples have no sequencing issues and a "
            "consistent protocol was followed for the sequencing of all samples. If this is "
            "not the case, deactivate this feature by activating this flag")
        ->group("Parameter Estimation");

    compare_subcmd
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Filtering");

    compare_subcmd
        ->add_flag(
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

    compare_subcmd
        ->add_flag(
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

    compare_subcmd
        ->add_flag(
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

    description = "Add extra step to carefully genotype sites.";
    auto* gt_opt = compare_subcmd->add_flag("--genotype", opt->genotype, description)
                       ->group("Consensus/Variant Calling");

    description = "(Intended for developers) Use coverage-oriented (local) genotyping "
                  "instead of the default ML path-oriented (global) approach.";
    compare_subcmd->add_flag("--local", opt->local_genotype, description)
        ->needs(gt_opt)
        ->group("Genotyping");

    description
        = "Minimum size of a cluster of hits between a read and a loci to consider "
          "the loci present";
    compare_subcmd
        ->add_option("-c,--min-cluster-size", opt->min_cluster_size, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Mapping");

    description = "Maximum number of kmers to average over when selecting the maximum "
                  "likelihood path";
    compare_subcmd->add_option("--kmer-avg", opt->max_num_kmers_to_avg, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Consensus/Variant Calling");

    description
        = "Hard threshold for the minimum allele coverage allowed when genotyping";
    compare_subcmd->add_option("-a", opt->min_allele_covg_gt, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

    description = "The minimum required total coverage for a site when genotyping";
    compare_subcmd->add_option("-s", opt->min_total_covg_gt, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

    description = "Minimum difference in coverage on a site required between the first "
                  "and second maximum likelihood path";
    compare_subcmd->add_option("-D", opt->min_diff_covg_gt, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Genotyping");

    description = "Minimum allele coverage, as a fraction of the expected coverage, "
                  "allowed when genotyping";
    compare_subcmd->add_option("-F", opt->min_allele_fraction_covg_gt, description)
        ->capture_default_str()
        ->type_name("INT")
        ->group("Genotyping");

    description = "When genotyping, assume that coverage on alternative alleles arises "
                  "as a result of an error process with rate -E.";
    compare_subcmd
        ->add_option("-E,--gt-error-rate", opt->genotyping_error_rate, description)
        ->capture_default_str()
        ->group("Genotyping");

    description = "Minimum genotype confidence (GT_CONF) required to make a call";
    compare_subcmd->add_option("-G,--gt-conf", opt->confidence_threshold, description)
        ->type_name("INT")
        ->capture_default_str()
        ->group("Genotyping");

    compare_subcmd
        ->add_flag("-K,--debugging-files", opt->keep_extra_debugging_files,
            "Keep extra debugging files. Warning: this might "
            "create thousands of files.")
        ->group("Debugging");

    compare_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    compare_subcmd->callback([opt]() { pandora_compare(*opt); });
}

int pandora_compare(CompareOptions& opt)
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
    // ==========
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

    auto samples = load_read_index(opt.reads_idx_file);
    std::vector<std::string> sample_names;
    for (const auto& sample_pair : samples) {
        sample_names.push_back(sample_pair.first);
    }

    auto pangraph = std::make_shared<pangenome::Graph>(sample_names);

    // map reads and get the sample pangraph for each sample
    for (uint32_t sample_id = 0; sample_id < samples.size(); ++sample_id) {
        const auto& sample = samples[sample_id];
        auto pangraph_sample = std::make_shared<pangenome::Graph>();

        const auto& sample_name = sample.first;
        const auto& sample_fpath = sample.second;

        // make output dir for this sample
        const auto sample_outdir { opt.outdir / sample_name };
        fs::create_directories(sample_outdir);

        BOOST_LOG_TRIVIAL(info) << "Constructing pangenome::Graph from read file "
                                << sample_fpath << " (this will take a while)";
        uint32_t covg = pangraph_from_read_file(sample, pangraph_sample, index,
            opt.max_diff, opt.error_rate, sample_outdir,
            opt.min_cluster_size, opt.genome_size, opt.illumina, opt.clean,
            opt.max_covg, opt.threads, opt.keep_extra_debugging_files);

        const auto pangraph_gfa { sample_outdir / "pandora.pangraph.gfa" };
        BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file " << pangraph_gfa;
        write_pangraph_gfa(pangraph_gfa, pangraph_sample);

        if (pangraph_sample->nodes.empty()) {
            BOOST_LOG_TRIVIAL(warning)
                << "Found no LocalPRGs in the reads for sample " << sample_name;
        }

        BOOST_LOG_TRIVIAL(info) << "Update LocalPRGs with hits";
        pangraph_sample->add_hits_to_kmergraphs(0);

        BOOST_LOG_TRIVIAL(info) << "Estimate parameters for kmer graph model";
        auto exp_depth_covg = estimate_parameters(pangraph_sample, sample_outdir,
            index.get_kmer_size(), opt.error_rate, covg, opt.binomial, 0,
            opt.do_not_auto_update_params);

        genotyping_options.add_exp_depth_covg(exp_depth_covg);

        if (genotyping_options.get_min_kmer_covg() == 0) {
            genotyping_options.set_min_kmer_covg(exp_depth_covg / 10);
        }

        BOOST_LOG_TRIVIAL(info) << "Find max likelihood PRG paths";
        auto sample_pangraph_size = pangraph_sample->nodes.size();
        Fastaq consensus_fq(true, true);
        for (auto c = pangraph_sample->nodes.begin();
             c != pangraph_sample->nodes.end();) {
            const auto& local_prg = index.get_prg_given_id(c->second->prg_id);
            vector<KmerNodePtr> kmp;
            vector<LocalNodePtr> lmp;
            local_prg->add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp,
                index.get_window_size(), opt.binomial, covg, opt.max_num_kmers_to_avg, 0,
                opt.min_absolute_gene_coverage, opt.min_relative_gene_coverage,
                opt.max_relative_gene_coverage);

            if (kmp.empty()) {
                c = pangraph_sample->remove_node(c->second);
                continue;
            }

            pangraph->add_node(local_prg);
            pangraph->add_hits_between_PRG_and_sample(
                c->second->prg_id, sample_name, kmp);

            ++c;
        }

        pangraph->copy_coverages_to_kmergraphs(*pangraph_sample, sample_id);

        consensus_fq.save(sample_outdir / "pandora.consensus.fq.gz");
        consensus_fq.clear();
        if (pangraph_sample->nodes.empty() and sample_pangraph_size > 0) {
            BOOST_LOG_TRIVIAL(warning)
                << "All LocalPRGs found were removed for sample " << sample_name
                << ". Is your genome_size accurate? Genome size is assumed to be "
                << opt.genome_size << " and can be updated with --genome_size";
        }

        // Note: pangraph_sample is destroyed here and as well as all Read information
        // (pangenome::Graph::reads) about the sample pangraph does not keep the read
        // information This is important since this is the heaviest information to keep
        // in compare pangraph has just coverage information and the consensus path for
        // each sample and PRG
    }

    // for each pannode in graph, find a best reference
    // and output a vcf and aligned fasta of sample paths through it
    BOOST_LOG_TRIVIAL(info) << "Multi-sample pangraph has " << pangraph->nodes.size()
                            << " nodes";

    // parallel region!

    // load vcf refs
    VCFRefs vcf_refs; // no need to control this variable - read only
    {
        if (!opt.vcf_refs_file.empty()) {
            vcf_refs.reserve(index.get_number_of_prgs());
            load_vcf_refs_file(opt.vcf_refs_file, vcf_refs);
        }
    }

    // shared variable - controlled by critical(vcf_ref_fa)
    Fastaq vcf_ref_fa(true, false);

    // shared variable - controlled by critical(VCFPathsToBeConcatenated)
    std::vector<fs::path> VCFPathsToBeConcatenated; // TODO: we can allocate the correct
                                                    // size of the vectors already here

    // shared variable - controlled by critical(VCFGenotypedPathsToBeConcatenated)
    std::vector<fs::path>
        VCFGenotypedPathsToBeConcatenated; // TODO: we can allocate the correct size of
                                           // the vectors already here

    // create the dir that will contain all vcfs
    const int nb_vcfs_per_dir = 4000;
    const auto vcfs_dir { opt.outdir / "VCFs" };
    fs::create_directories(vcfs_dir);
    // create the dirs for the VCFs
    for (uint32_t i = 0; i <= pangraph->nodes.size() / nb_vcfs_per_dir; ++i) {
        fs::create_directories(vcfs_dir / int_to_string(i + 1));
    }

    // create the dirs for the VCFs genotyped, if genotyping should be done
    const auto vcfs_genotyped_dirs { opt.outdir / "VCFs_genotyped" };
    if (opt.genotype) {
        for (uint32_t i = 0; i <= pangraph->nodes.size() / nb_vcfs_per_dir; ++i) {
            fs::create_directories(vcfs_genotyped_dirs / int_to_string(i + 1));
        }
    }

    // transforms to a vector to parallelize this
    // TODO: use OMP task instead?
    std::vector<std::pair<NodeId, std::shared_ptr<pangenome::Node>>>
        pangraphNodesAsVector;
    pangraphNodesAsVector.reserve(pangraph->nodes.size());
    for (auto pan_id_to_node_mapping = pangraph->nodes.begin();
         pan_id_to_node_mapping != pangraph->nodes.end(); ++pan_id_to_node_mapping) {
        // todo: use emplace_back?
        pangraphNodesAsVector.push_back(*pan_id_to_node_mapping);
    }

#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 1) default(shared)
    for (uint32_t pangraph_node_index = 0;
         pangraph_node_index < pangraphNodesAsVector.size(); ++pangraph_node_index) {
        const auto& pangraph_node_entry = pangraphNodesAsVector[pangraph_node_index];
        pangenome::Node& pangraph_node = *pangraph_node_entry.second;

        const auto& prg_id = pangraph_node.prg_id;

        const bool valid_prg_id = index.get_number_of_prgs() >= prg_id;
        if (!valid_prg_id) {
            fatal_error("Error reading PanRG: a PRG has an invalid ID (", prg_id,
                "), >= than the number of PRGs (", index.get_number_of_prgs(), ") in the PanRG");
        }
        const auto& prg_ptr = index.get_prg_given_id(prg_id);

        const auto vcf_reference_path
            = pangraph->infer_node_vcf_reference_path(pangraph_node, prg_ptr,
                index.get_window_size(), vcf_refs, opt.max_num_kmers_to_avg);

#pragma omp critical(vcf_ref_fa)
        {
            vcf_ref_fa.add_entry(
                prg_ptr->name, prg_ptr->string_along_path(vcf_reference_path), "");
        }

        // output the vcf for this sample
        VCF vcf(&genotyping_options);

        // add all samples to the vcf
        vcf.add_samples(sample_names);

        // build the vcf
        pangraph_node.construct_multisample_vcf(
            vcf, vcf_reference_path, prg_ptr, index.get_window_size());

        // save the vcf to disk
        uint32_t dir = pangraph_node_index / nb_vcfs_per_dir
            + 1; // get the good dir for this sample vcf
        const auto vcf_path { vcfs_dir / int_to_string(dir)
            / (prg_ptr->name + ".vcf") };
        vcf.save(vcf_path, true, false);

// add the vcf path to VCFPathsToBeConcatenated to concatenate after
#pragma omp critical(VCFPathsToBeConcatenated)
        {
            VCFPathsToBeConcatenated.push_back(vcf_path);
        }

        if (opt.genotype) {
            BOOST_LOG_TRIVIAL(info) << "Genotyping VCF...";
            vcf.genotype(opt.local_genotype);
            BOOST_LOG_TRIVIAL(info) << "Finished genotyping VCF";

            // save the genotyped vcf to disk
            const auto vcf_genotyped_path { vcfs_genotyped_dirs / int_to_string(dir)
                / (prg_ptr->name + "_genotyped.vcf") };
            vcf.save(vcf_genotyped_path, false, true);

// add the genotyped vcf path to VCFGenotypedPathsToBeConcatenated to concatenate after
#pragma omp critical(VCFGenotypedPathsToBeConcatenated)
            {
                VCFGenotypedPathsToBeConcatenated.push_back(vcf_genotyped_path);
            }
        }
    }

    // generate all the multisample files
    vcf_ref_fa.save(opt.outdir / "pandora_multisample.vcf_ref.fa");
    VCF::concatenate_VCFs(
        VCFPathsToBeConcatenated, opt.outdir / "pandora_multisample_consensus.vcf");
    if (opt.genotype) {
        VCF::concatenate_VCFs(VCFGenotypedPathsToBeConcatenated,
            opt.outdir / "pandora_multisample_genotyped.vcf");
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    BOOST_LOG_TRIVIAL(info) << "Output matrix";
    pangraph->save_matrix(opt.outdir / "pandora_multisample.matrix", sample_names);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "No LocalPRGs found to compare samples on. "
            << "Is your genome_size accurate? Genome size is assumed to be "
            << opt.genome_size << " and can be updated with --genome-size";
    }

    BOOST_LOG_TRIVIAL(info) << "Done!";
    return 0;
}
