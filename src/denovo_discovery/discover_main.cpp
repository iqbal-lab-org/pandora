//
// Created by michael on 26/9/20.
//

#include "denovo_discovery/discover_main.h"
#include "denovo_discovery/racon.h"
#include "denovo_discovery/denovo_record.h"
#include "cli_helpers.h"

void setup_discover_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<DiscoverOptions>();
    auto* discover_subcmd = app.add_subcommand("discover",
        "Quasi-map reads from multiple samples to a pandora index, infer the "
        "sequence of present loci in the sample and discover novel variants.");

    discover_subcmd
        ->add_option("<TARGET>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
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

    description
        = "Minimum proportion of overlap between two clusters of hits to consider "
          "them conflicting. Only one of the conflicting clusters will be kept.";
    discover_subcmd
        ->add_option("--min-cluster-overlap", opt->conflicting_clusters_overlap_threshold, description)
        ->capture_default_str()
        ->type_name("FLOAT")
        ->group("Mapping");

    description
        = "When two clusters of hits are conflicting, the one with highest number of unique minimisers "
          "will be kept. However, if the difference between the number of unique minimisers is too small, "
          "less than this parameter, then we will prefer the cluster that has higher target coverage.";
    discover_subcmd
        ->add_option("--cluster-mini-tolerance", opt->conflicting_clusters_minimiser_tolerance, description)
        ->capture_default_str()
        ->type_name("FLOAT")
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
        ->add_flag("--bin", opt->binomial,
            "Use binomial model for kmer coverages [default: negative binomial]")
        ->group("Parameter Estimation");

    discover_subcmd
        ->add_flag("--dont-auto-update-params", opt->do_not_auto_update_params,
            "By default, pandora automatically updates error rate and kmer coverage model parameters based "
            "on the mapping of the previous sample. This could potentially generate "
            "more accurate results if your samples have no sequencing issues and a "
            "consistent protocol was followed for the sequencing of all samples. If this is "
            "not the case, deactivate this feature by activating this flag")
        ->group("Parameter Estimation");

    discover_subcmd
        ->add_option("--max-covg", opt->max_covg, "Maximum coverage of reads to accept")
        ->capture_default_str()
        ->type_name("INT")
        ->group("Filtering");

    discover_subcmd
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

    discover_subcmd
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

    discover_subcmd
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

    discover_subcmd
        ->add_flag("-K,--debugging-files", opt->keep_extra_debugging_files,
            "Keep extra debugging files. Warning: this might "
            "create thousands of files.")
        ->group("Debugging");

    discover_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    discover_subcmd
        ->add_option(
            "-r,--rng-seed", opt->rng_seed, "RNG seed, an int>0 to force deterministic "
                                            "mapping when multiple optimal mappings are "
                                            "possible. To be avoided except in "
                                            "debugging/investigation scenarios. A value "
                                            "of 0 will be interpreted as no seed given "
                                            "and mapping will not be deterministic.")
        ->capture_default_str()
        ->group("Debugging");

    // Set the function that will be called when this subcommand is issued.
    discover_subcmd->callback([opt]() { pandora_discover(*opt); });
}

void pandora_discover_core(const SampleData& sample, Index &index, DiscoverOptions& opt)
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
        = pangraph_from_read_file(sample, pangraph, index, opt.max_diff, opt.error_rate, sample_outdir,
        opt.min_cluster_size, opt.genome_size, opt.max_covg,
        opt.conflicting_clusters_overlap_threshold, opt.conflicting_clusters_minimiser_tolerance,
        opt.threads, opt.keep_extra_debugging_files, opt.rng_seed);

    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(warning)
            << "[Sample " << sample_name << "] "
            << "Found no LocalPRGs in the reads for sample " << sample_name;
    }

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Updating error rate...";
    estimate_parameters(pangraph, sample_outdir, index.get_kmer_size(), opt.error_rate,
        covg, opt.binomial, 0, opt.do_not_auto_update_params);

    BOOST_LOG_TRIVIAL(info) << "[Sample " << sample_name << "] "
                            << "Find PRG paths and discover novel alleles...";

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

    // TODO: this can be memory heavy, fix
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
        const auto &prg = index.get_prg_given_id(pangraph_node->prg_id);
        prg->add_consensus_path_to_fastaq(consensus_fq,
             pangraph_node, kmp, lmp, index.get_window_size(), opt.binomial, covg,
             opt.max_num_kmers_to_avg, 0,
             opt.min_absolute_gene_coverage, opt.min_relative_gene_coverage,
             opt.max_relative_gene_coverage);

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
                kmer_graph_dir / (pangraph_node->get_name() + ".kg.gfa"), prg);
        }

        // builds a mem_fd with the locus reads
        const std::pair<int, std::string> locus_reads_fd_and_filepath = build_memfd(locus_to_reads[locus]);

        if (opt.keep_extra_debugging_files) {
            // build the reads file on disk
            build_file((denovo_outdir / (locus + ".reads.fa")).string(),
                locus_to_reads[locus]);
        }

        const std::string lmp_seq = prg->string_along_path(lmp);
        Racon racon(opt.illumina, index.get_kmer_size(), locus, lmp_seq,
            denovo_outdir, locus_reads_fd_and_filepath.second, 10, opt.keep_extra_debugging_files);
        const std::string &polished_sequence = racon.get_polished_sequence();

        close(locus_reads_fd_and_filepath.first);

        std::string denovo_sequence;
        denovo_sequence.reserve(polished_sequence.size() + 1024);
        denovo_sequence
                .append(">")
                .append(locus)
                .append(" sample=")
                .append(sample_name)
                .append(" denovo_sequence\n")
                .append(polished_sequence);

#pragma omp critical(all_denovo_sequences)
        {
            all_denovo_sequences.push_back(denovo_sequence);
        }

        const std::vector<DenovoVariantRecord> denovo_variants =
            DenovoVariantRecord::get_variants_from_pair_of_sequences(lmp_seq, polished_sequence);

        const bool we_have_denovo_variants = denovo_variants.size() > 0;
        if (we_have_denovo_variants) {
            std::string denovo_paths_description;
            {
                std::stringstream denovo_paths_out_core_file_ss;
                denovo_paths_out_core_file_ss << locus << std::endl;
                const std::string lmp_as_string = LocalNode::to_string_vector(lmp);
                denovo_paths_out_core_file_ss << lmp_as_string << std::endl;
                denovo_paths_out_core_file_ss << denovo_variants.size()
                                              << " denovo variants for this locus" << std::endl;
                for (const DenovoVariantRecord &denovo_variant : denovo_variants) {
                    denovo_paths_out_core_file_ss << denovo_variant.to_string() << std::endl;
                }

                denovo_paths_description = denovo_paths_out_core_file_ss.str();
            }

#pragma omp critical(denovo_paths_out_core_file)
            {
                ++number_of_loci_with_denovo_variants;
                denovo_paths_out_core_file << denovo_paths_description;
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

    if (!opt.keep_extra_debugging_files) {
        fs::remove_all(denovo_outdir);
    }

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

    if (opt.illumina) {
        opt.error_rate = 0.001;
    }

    BOOST_LOG_TRIVIAL(info) << "Loading Index...";
    Index index = Index::load(opt.index_file.string());
    BOOST_LOG_TRIVIAL(info) << "Index loaded successfully!";

    if (opt.illumina and opt.max_diff > 200) {
        opt.max_diff = 2 * index.get_kmer_size() + 1;
    }

    auto samples = load_read_index(opt.reads_idx_file);

    // for each sample, run pandora discover
    for (const SampleData& sample : samples) {
        pandora_discover_core(sample, index, opt);
    }

    concatenate_all_denovo_files(samples, opt.outdir);
    BOOST_LOG_TRIVIAL(info) << "All done!";
    return 0;
}
