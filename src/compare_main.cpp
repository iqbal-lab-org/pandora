#include "compare_main.h"

void setup_compare_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<CompareOptions>();
    auto* compare_subcmd = app.add_subcommand("compare",
        "Quasi-map reads from multiple samples to an indexed PRG, infer the "
        "sequence of present loci in each sample, and call variants between the "
        "samples.");

    compare_subcmd
        ->add_option("<TARGET>", opt->prgfile, "An indexed PRG file (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    std::string description
        = "A tab-delimited file where each line is a sample identifier followed by "
          "the path to the fast{a,q} of reads for that sample";
    compare_subcmd->add_option("<QUERY_IDX>", opt->reads_idx_file, description)
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    compare_subcmd
        ->add_option(
            "-w", opt->window_size, "Window size for (w,k)-minimizers (must be <=k)")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    compare_subcmd->add_option("-k", opt->kmer_size, "K-mer size for (w,k)-minimizers")
        ->type_name("INT")
        ->capture_default_str()
        ->group("Indexing");

    compare_subcmd
        ->add_option("-o,--outdir", opt->outdir, "Directory to write output files to")
        ->type_name("DIR")
        ->capture_default_str()
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
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
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

    compare_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    compare_subcmd->callback([opt]() { pandora_compare(*opt); });
}

std::vector<std::pair<SampleIdText, SampleFpath>> load_read_index(
    const std::string& read_index_fpath)
{
    std::map<SampleIdText, SampleFpath> samples;
    std::string name, reads_path, line;
    std::ifstream myfile(read_index_fpath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            std::istringstream linestream(line);
            if (std::getline(linestream, name, '\t')) {
                linestream >> reads_path;
                if (samples.find(name) != samples.end()) {
                    BOOST_LOG_TRIVIAL(warning)
                        << "Warning: non-unique sample ids given! Only the last "
                           "of these will be kept";
                }
                samples[name] = reads_path;
            }
        }
    } else {
        BOOST_LOG_TRIVIAL(error)
            << "Unable to open read index file " << read_index_fpath;
        exit(1);
    }
    BOOST_LOG_TRIVIAL(info) << "Finished loading " << samples.size()
                            << " samples from read index";
    return std::vector<std::pair<SampleIdText, SampleFpath>>(
        samples.begin(), samples.end());
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

    BOOST_LOG_TRIVIAL(info) << "Loading Index and LocalPRGs from file...";
    auto index = std::make_shared<Index>();
    index->load(opt.prgfile, opt.window_size, opt.kmer_size);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);
    load_PRG_kmergraphs(prgs, opt.window_size, opt.kmer_size, opt.prgfile);

    BOOST_LOG_TRIVIAL(info) << "Loading read index file...";
    auto samples = load_read_index(opt.reads_idx_file);
    std::vector<std::string> sample_names;
    for (const auto& sample_pair : samples) {
        sample_names.push_back(sample_pair.first);
    }

    auto pangraph = std::make_shared<pangenome::Graph>(sample_names);

    // for each sample, run pandora to get the sample pangraph
    for (uint32_t sample_id = 0; sample_id < samples.size(); ++sample_id) {
        const auto& sample = samples[sample_id];
        auto pangraph_sample = std::make_shared<pangenome::Graph>();

        const auto& sample_name = sample.first;
        const auto& sample_fpath = sample.second;

        // make output dir for this sample
        auto sample_outdir = opt.outdir;
        sample_outdir.append("/").append(sample_name);
        fs::create_directories(sample_outdir);

        BOOST_LOG_TRIVIAL(info) << "Constructing pangenome::Graph from read file "
                                << sample_fpath << " (this will take a while)";
        uint32_t covg = pangraph_from_read_file(sample_fpath, pangraph_sample, index,
            prgs, opt.window_size, opt.kmer_size, opt.max_diff, opt.error_rate,
            opt.min_cluster_size, opt.genome_size, opt.illumina, opt.clean,
            opt.max_covg, opt.threads);
        BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file " << sample_outdir
                                << "/pandora.pangraph.gfa";
        write_pangraph_gfa(sample_outdir + "/pandora.pangraph.gfa", pangraph_sample);

        if (pangraph_sample->nodes.empty()) {
            BOOST_LOG_TRIVIAL(warning)
                << "Found no LocalPRGs in the reads for sample " << sample_name;
        }

        BOOST_LOG_TRIVIAL(info) << "Update LocalPRGs with hits";
        pangraph_sample->add_hits_to_kmergraphs(prgs, 0);

        BOOST_LOG_TRIVIAL(info) << "Estimate parameters for kmer graph model";
        auto exp_depth_covg = estimate_parameters(pangraph_sample, sample_outdir,
            opt.kmer_size, opt.error_rate, covg, opt.binomial, 0);
        genotyping_options.add_exp_depth_covg(exp_depth_covg);

        if (genotyping_options.get_min_kmer_covg() == 0) {
            genotyping_options.set_min_kmer_covg(exp_depth_covg / 10);
        }

        BOOST_LOG_TRIVIAL(info) << "Find max likelihood PRG paths";
        auto sample_pangraph_size = pangraph_sample->nodes.size();
        Fastaq consensus_fq(true, true);
        for (auto c = pangraph_sample->nodes.begin();
             c != pangraph_sample->nodes.end();) {
            const LocalPRG& local_prg = *prgs[c->second->prg_id];
            vector<KmerNodePtr> kmp;
            vector<LocalNodePtr> lmp;
            local_prg.add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp,
                opt.window_size, opt.binomial, covg, opt.max_num_kmers_to_avg, 0);

            if (kmp.empty()) {
                c = pangraph_sample->remove_node(c->second);
                continue;
            }

            pangraph->add_node(prgs[c->second->prg_id]);
            pangraph->add_hits_between_PRG_and_sample(
                c->second->prg_id, sample_name, kmp);

            ++c;
        }

        pangraph->copy_coverages_to_kmergraphs(*pangraph_sample, sample_id);

        consensus_fq.save(sample_outdir + "/pandora.consensus.fq.gz");
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
            vcf_refs.reserve(prgs.size());
            load_vcf_refs_file(opt.vcf_refs_file, vcf_refs);
        }
    }

    // shared variable - controlled by critical(vcf_ref_fa)
    Fastaq vcf_ref_fa(true, false);

    // shared variable - controlled by critical(VCFPathsToBeConcatenated)
    std::vector<std::string>
        VCFPathsToBeConcatenated; // TODO: we can allocate the correct size of the
                                  // vectors already here

    // shared variable - controlled by critical(VCFGenotypedPathsToBeConcatenated)
    std::vector<std::string>
        VCFGenotypedPathsToBeConcatenated; // TODO: we can allocate the correct size of
                                           // the vectors already here

    // create the dir that will contain all vcfs
    const int nb_vcfs_per_dir = 4000;
    auto vcfs_dir = opt.outdir + "/VCFs";
    fs::create_directories(vcfs_dir);
    // create the dirs for the VCFs
    for (uint32_t i = 0; i <= pangraph->nodes.size() / nb_vcfs_per_dir; ++i) {
        fs::create_directories(vcfs_dir + "/" + int_to_string(i + 1));
    }

    // create the dirs for the VCFs genotyped, if genotyping should be done
    auto vcfs_genotyped_dirs = opt.outdir + "/VCFs_genotyped";
    if (opt.genotype) {
        for (uint32_t i = 0; i <= pangraph->nodes.size() / nb_vcfs_per_dir; ++i) {
            fs::create_directories(vcfs_genotyped_dirs + "/" + int_to_string(i + 1));
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

#pragma omp parallel for num_threads(opt.threads) schedule(dynamic, 1)
    for (uint32_t pangraph_node_index = 0;
         pangraph_node_index < pangraphNodesAsVector.size(); ++pangraph_node_index) {
        const auto& pangraph_node_entry = pangraphNodesAsVector[pangraph_node_index];
        pangenome::Node& pangraph_node = *pangraph_node_entry.second;

        const auto& prg_id = pangraph_node.prg_id;
        assert(prgs.size() > prg_id);
        const auto& prg_ptr = prgs[prg_id];

        const auto vcf_reference_path
            = pangraph->infer_node_vcf_reference_path(pangraph_node, prg_ptr,
                opt.window_size, vcf_refs, opt.max_num_kmers_to_avg);

#pragma omp critical(vcf_ref_fa)
        {
            vcf_ref_fa.add_entry(
                prg_ptr->name, prg_ptr->string_along_path(vcf_reference_path), "");
        }

        // output the vcf for this sample
        VCF vcf(&genotyping_options);

        // add all samples to the vcf
        vcf.add_samples(sample_names);
        assert(vcf.samples.size() == samples.size());

        // build the vcf
        pangraph_node.construct_multisample_vcf(
            vcf, vcf_reference_path, prg_ptr, opt.window_size);

        // save the vcf to disk
        uint32_t dir = pangraph_node_index / nb_vcfs_per_dir
            + 1; // get the good dir for this sample vcf
        auto vcf_path
            = vcfs_dir + "/" + int_to_string(dir) + "/" + prg_ptr->name + ".vcf";
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
            auto vcf_genotyped_path = vcfs_genotyped_dirs + "/" + int_to_string(dir)
                + "/" + prg_ptr->name + "_genotyped.vcf";
            vcf.save(vcf_genotyped_path, false, true);

// add the genotyped vcf path to VCFGenotypedPathsToBeConcatenated to concatenate after
#pragma omp critical(VCFGenotypedPathsToBeConcatenated)
            {
                VCFGenotypedPathsToBeConcatenated.push_back(vcf_genotyped_path);
            }
        }
    }

    // generate all the multisample files
    vcf_ref_fa.save(opt.outdir + "/pandora_multisample.vcf_ref.fa");
    VCF::concatenate_VCFs(
        VCFPathsToBeConcatenated, opt.outdir + "/pandora_multisample_consensus.vcf");
    if (opt.genotype) {
        VCF::concatenate_VCFs(VCFGenotypedPathsToBeConcatenated,
            opt.outdir + "/pandora_multisample_genotyped.vcf");
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    BOOST_LOG_TRIVIAL(info) << "Output matrix";
    pangraph->save_matrix(opt.outdir + "/pandora_multisample.matrix", sample_names);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(error)
            << "No LocalPRGs found to compare samples on. "
            << "Is your genome_size accurate? Genome size is assumed to be "
            << opt.genome_size << " and can be updated with --genome_size";
    }

    // clear up
    index->clear();

    // current date/time based on current system
    BOOST_LOG_TRIVIAL(info) << "Done!";
    return 0;
}
