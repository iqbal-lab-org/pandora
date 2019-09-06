/*
 * C++ Program to call variation between a set of samples. 
 * To do this, runs pandora map process on each sample/readfile to find the best gene sequence for 
 * each gene inferred present in that sample. Then compares the paths found through each gene and 
 * outputs a vcf and aligned fasta for each one. 
 */
/*
 * QUESTIONS:
 * How do I handle multiple occurrences of a gene in a sample? I would like to output all the versions called against my new reference
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <tuple>
#include <functional>
#include <cctype>
#include <fstream>
#include <algorithm>
#include <map>
#include <cassert>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "noise_filtering.h"
#include "estimate_parameters.h"


using std::set;
using std::vector;


static void show_compare_usage() {
    std::cerr << "Usage: pandora compare -p PRG_FILE -r READ_INDEX -o OUTDIR <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-p,--prg_file PRG_FILE\t\tSpecify a fasta-style prg file\n"
              << "\t-r,--read_index READ_INDEX\tSpecify a file with a line per sample\n"
              << "\t\t\t\t\tsample_id <tab> filepath to reads in fasta/q format\n"
              << "\t-o,--outdir OUTDIR\tSpecify directory of output\n"
              << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, default 14\n"
              << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
              << "\t-m,--max_diff INT\t\tMaximum distance between consecutive hits within a cluster, default 250 (bps)\n"
              << "\t-e,--error_rate FLOAT\t\tEstimated error rate for reads, default 0.11\n"
              << "\t-t T\t\t\t\tNumber of threads, default 1\n"
              << "\t--genome_size\tNUM_BP\tEstimated length of genome, used for coverage estimation\n"
              << "\t--vcf_refs REF_FASTA\t\tA fasta file with an entry for each LocalPRG giving reference sequence for\n"
              << "\t\t\t\t\tVCF. Must have a perfect match in the graph and the same name as the graph\n"
              << "\t--illumina\t\t\tData is from illumina rather than nanopore, so is shorter with low error rate\n"
              << "\t--clean\t\t\tAdd a step to clean and detangle the pangraph\n"
              << "\t--bin\t\t\tUse binomial model for kmer coverages, default is negative binomial\n"
              << "\t--max_covg\t\t\tMaximum average coverage from reads to accept\n"
              << "\t--genotype\t\t\tAdd extra step to carefully genotype sites\n"
              << "\t--log_level\t\t\tdebug,[info],warning,error\n"
              << std::endl;
}

using SampleIdText = std::string;
using SampleFpath = std::string;

std::vector<std::pair<SampleIdText, SampleFpath>> load_read_index(const std::string &read_index_fpath) {
    std::map<SampleIdText, SampleFpath> samples;
    std::string name, reads_path, line;
    std::ifstream myfile(read_index_fpath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            std::istringstream linestream(line);
            if (std::getline(linestream, name, '\t')) {
                linestream >> reads_path;
                if (samples.find(name) != samples.end()) {
                    std::cout << "Warning: non-unique sample ids given! Only the last of these will be kept"
                              << std::endl;
                }
                samples[name] = reads_path;
            }
        }
    } else {
        std::cerr << "Unable to open read index file " << read_index_fpath << std::endl;
        exit(1);
    }
    std::cout << now() << "Finished loading " << samples.size() << " samples from read index" << std::endl;
    return std::vector<std::pair<SampleIdText, SampleFpath>>(samples.begin(), samples.end());
}

int pandora_compare(int argc, char *argv[]) {
    // if not enough arguments, print usage
    if (argc < 7) {
        show_compare_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    std::string prgfile, read_index_fpath, outdir = "pandora", vcf_refs_file, log_level="info";
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000, max_covg = 300, max_num_kmers_to_average=100,
            min_allele_covg_gt = 0, min_total_covg_gt = 0, min_diff_covg_gt = 0, min_kmer_covg=0, threads=1; // default parameters
    uint16_t confidence_threshold = 1;
    int max_diff = 250;
    float e_rate = 0.11, min_allele_fraction_covg_gt = 0, genotyping_error_rate=0.01;
    bool illumina = false, clean = false, bin = false, genotype = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_compare_usage();
            return 0;
        } else if ((arg == "-p") || (arg == "--prg_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prgfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--prg_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-r") || (arg == "--read_index")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                read_index_fpath = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--read_index option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-o") || (arg == "--outdir")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                outdir = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--outdir option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-w") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                w = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-w option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-k") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                k = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-k option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-m") || (arg == "--max_diff")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_diff = strtol(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_diff option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-c") || (arg == "--min_cluster_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_cluster_size = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_cluster_size option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-e") || (arg == "--error_rate")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                e_rate = strtof(argv[++i], nullptr); // Increment 'i' so we don't get the argument as the next argv[i].
                if (e_rate < 0.01) {
                    illumina = true;
                }
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--error_rate option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genome_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                genome_size = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--genome_size option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "--vcf_refs") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                vcf_refs_file = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--vcf_refs option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--illumina")) {
            illumina = true;
        } else if ((arg == "--clean")) {
            clean = true;
        } else if ((arg == "--bin")) {
            bin = true;
        } else if ((arg == "--max_covg")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_covg = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_covg option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--max_num_kmers_to_average")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_num_kmers_to_average = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_num_kmers_to_average option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_allele_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_allele_covg_gt = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_allele_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_total_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_total_covg_gt = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_total_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_diff_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_diff_covg_gt = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_diff_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_allele_fraction_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_allele_fraction_covg_gt = strtof(argv[++i], nullptr); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_allele_fraction_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotyping_error_rate")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                genotyping_error_rate = strtof(argv[++i], nullptr); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--genotyping_error_rate option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--confidence_threshold")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                confidence_threshold = (uint16_t)strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--confidence_threshold option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotype")) {
            genotype = true;
        } else if ((arg == "--log_level")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                log_level = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--log_level option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-t") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                threads = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-t option requires one argument." << std::endl;
                return 1;
            }
        }
        else {
            std::cerr << argv[i] << " could not be attributed to any parameter" << std::endl;
        }
    }

    assert(w <= k);
    assert(not prgfile.empty());
    if (illumina and e_rate > 0.1) {
        e_rate = 0.001;
    }
    if (illumina and max_diff > 200) {
        max_diff = 2*k + 1;
    }

    //then run the programme...
    std::cout << "START: " << now() << std::endl;
    std::cout << "\nUsing parameters: " << std::endl;
    std::cout << "\tprgfile\t\t" << prgfile << std::endl;
    std::cout << "\tread_index_fpath\t" << read_index_fpath << std::endl;
    std::cout << "\toutdir\t" << outdir << std::endl;
    std::cout << "\tw\t\t" << w << std::endl;
    std::cout << "\tk\t\t" << k << std::endl;
    std::cout << "\tmax_diff\t" << max_diff << std::endl;
    std::cout << "\terror_rate\t" << e_rate << std::endl;
    std::cout << "\tthreads\t" << threads << std::endl;
    std::cout << "\tvcf_refs\t" << vcf_refs_file << std::endl;
    std::cout << "\tillumina\t" << illumina << std::endl;
    std::cout << "\tclean\t" << clean << std::endl;
    std::cout << "\tbin\t" << bin << std::endl << std::endl;
    std::cout << "\tmax_covg\t" << max_covg << std::endl;
    std::cout << "\tgenotype\t" << genotype << std::endl;
    std::cout << "\tlog_level\t" << log_level << std::endl << std::endl;

    auto g_log_level{boost::log::trivial::info};
    if (log_level == "debug")
        g_log_level = boost::log::trivial::debug;
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= g_log_level);

    std::cout << now() << "Loading Index and LocalPRGs from file" << std::endl;
    auto index = std::make_shared<Index>();
    index->load(prgfile, w, k);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    // load read index
    std::cout << now() << "Loading read index file " << read_index_fpath << std::endl;
    auto samples = load_read_index(read_index_fpath);
    std::vector<std::string> sample_names;
    for (const auto &sample_pair : samples) sample_names.push_back(sample_pair.first);

    auto pangraph = std::make_shared<pangenome::Graph>(sample_names);
    std::vector<uint32_t> exp_depth_covgs;

    // for each sample, run pandora to get the sample pangraph
    for (uint32_t sample_id = 0; sample_id < samples.size(); ++sample_id) {
        const auto &sample = samples[sample_id];
        auto pangraph_sample = std::make_shared<pangenome::Graph>();

        const auto &sample_name = sample.first;
        const auto &sample_fpath = sample.second;

        // make output dir for this sample
        auto sample_outdir = outdir;
        sample_outdir.append("/").append(sample_name);
        fs::create_directories(sample_outdir);

        // construct the pangraph for this sample
        BOOST_LOG_TRIVIAL(info) << "Constructing pangenome::Graph from read file "
                                << sample_fpath
                                << " (this will take a while)";
        uint32_t covg = pangraph_from_read_file(sample_fpath,
                                                pangraph_sample,
                                                index, prgs, w, k,
                                                max_diff, e_rate,
                                                min_cluster_size, genome_size, illumina, clean, max_covg, threads);
        BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file "
                                << sample_outdir << "/pandora.pangraph.gfa";
        write_pangraph_gfa(sample_outdir + "/pandora.pangraph.gfa", pangraph_sample);

        if (pangraph_sample->nodes.empty()) {
            BOOST_LOG_TRIVIAL(warning) << "Found no LocalPRGs in the reads for sample " << sample_name;
        }

        BOOST_LOG_TRIVIAL(info) << "Update LocalPRGs with hits";
        pangraph_sample->add_hits_to_kmergraphs(prgs, 0);

        BOOST_LOG_TRIVIAL(info) << "Estimate parameters for kmer graph model";
        auto exp_depth_covg = estimate_parameters(pangraph_sample, sample_outdir, k, e_rate, covg, bin, 0);
        exp_depth_covgs.push_back(exp_depth_covg);

        if (min_kmer_covg == 0)
            min_kmer_covg = exp_depth_covg/10;

        BOOST_LOG_TRIVIAL(info) << "Find max likelihood PRG paths";
        auto sample_pangraph_size = pangraph_sample->nodes.size();
        Fastaq consensus_fq(true, true);
        for (auto c = pangraph_sample->nodes.begin(); c != pangraph_sample->nodes.end();) {
            const LocalPRG &local_prg = *prgs[c->second->prg_id];
            vector<KmerNodePtr> kmp;
            vector<LocalNodePtr> lmp;
            local_prg.add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp, w, bin, covg,
                                                   max_num_kmers_to_average, 0);

            if (kmp.empty()) {
                c = pangraph_sample->remove_node(c->second);
                continue;
            }

            pangraph->add_node(prgs[c->second->prg_id]);
            pangraph->add_hits_between_PRG_and_sample(c->second->prg_id, sample_name, kmp);

            ++c;
        }

        pangraph->copy_coverages_to_kmergraphs(*pangraph_sample, sample_id);

        consensus_fq.save(sample_outdir + "/pandora.consensus.fq.gz");
        consensus_fq.clear();
        if (pangraph_sample->nodes.empty() and sample_pangraph_size > 0) {
            std::cout << "WARNING: All LocalPRGs found were removed for sample " << sample_name
                      << ". Is your genome_size accurate? Genome size is assumed to be " << genome_size
                      << " and can be updated with --genome_size" << std::endl;
        }


        // Note: pangraph_sample is destroyed here and as well as all Read information (pangenome::Graph::reads) about the sample
        // pangraph does not keep the read information
        // This is important since this is the heaviest information to keep in compare
        // pangraph has just coverage information and the consensus path for each sample and PRG
    }

    // for each pannode in graph, find a best reference
    // and output a vcf and aligned fasta of sample paths through it
    BOOST_LOG_TRIVIAL(info) << "Multi-sample pangraph has "
                            << pangraph->nodes.size() << " nodes";







    //parallel region!

    // load vcf refs
    VCFRefs vcf_refs; //no need to control this variable - read only
    {
        std::string vcf_ref;
        if (!vcf_refs_file.empty()) {
            vcf_refs.reserve(prgs.size());
            load_vcf_refs_file(vcf_refs_file, vcf_refs);
        }
    }

    //shared variable - controlled by critical(vcf_ref_fa)
    Fastaq vcf_ref_fa(true, false);

    //shared variable - controlled by critical(VCFPathsToBeConcatenated)
    std::vector<std::string> VCFPathsToBeConcatenated; //TODO: we can allocate the correct size of the vectors already here

    //shared variable - controlled by critical(VCFGenotypedPathsToBeConcatenated)
    std::vector<std::string> VCFGenotypedPathsToBeConcatenated; //TODO: we can allocate the correct size of the vectors already here

    //create the dir that will contain all vcfs
    const int nbOfVCFsPerDir=4000;
    auto VCFsDirs = outdir + "/VCFs";
    fs::create_directories(VCFsDirs);
    //create the dirs for the VCFs
    for (uint32_t i = 0; i <= pangraph->nodes.size()/nbOfVCFsPerDir; ++i)
        fs::create_directories(VCFsDirs + "/" + int_to_string(i + 1));

    //create the dirs for the VCFs genotyped, if flag is passed
    auto VCFsGenotypedDirs = outdir + "/VCFs_genotyped";
    if (genotype) {
        for (uint32_t i = 0; i <= pangraph->nodes.size()/nbOfVCFsPerDir; ++i)
            fs::create_directories(VCFsGenotypedDirs + "/" + int_to_string(i + 1));
    }




    //transforms to a vector to parallelize this
    //TODO: use OMP task instead?
    std::vector<std::pair<NodeId, std::shared_ptr<pangenome::Node>>> pangraphNodesAsVector;
    pangraphNodesAsVector.reserve(pangraph->nodes.size());
    for (auto pan_id_to_node_mapping = pangraph->nodes.begin(); pan_id_to_node_mapping != pangraph->nodes.end(); ++pan_id_to_node_mapping)
        pangraphNodesAsVector.push_back(*pan_id_to_node_mapping);

    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint32_t pangraph_node_index=0; pangraph_node_index < pangraphNodesAsVector.size(); ++pangraph_node_index) {
        const auto &pangraph_node_entry = pangraphNodesAsVector[pangraph_node_index];
        BOOST_LOG_TRIVIAL(debug) << "Consider next node";
        pangenome::Node &pangraph_node = *pangraph_node_entry.second;

        const auto &prg_id = pangraph_node.prg_id;
        assert(prgs.size() > prg_id);
        const auto &prg_ptr = prgs[prg_id];

        const auto vcf_reference_path = pangraph->infer_node_vcf_reference_path(pangraph_node, prg_ptr, w, vcf_refs,
                                                                                max_num_kmers_to_average);

        #pragma omp critical(vcf_ref_fa)
        {
            vcf_ref_fa.add_entry(prg_ptr->name, prg_ptr->string_along_path(vcf_reference_path), "");
        }

        //output the vcf for this sample
        VCF vcf;

        //add all samples to the vcf
        vcf.add_samples(sample_names);
        assert( vcf.samples.size() == samples.size());

        //build the vcf
        pangraph_node.construct_multisample_vcf(vcf, vcf_reference_path, prg_ptr, w, min_kmer_covg);

        //save the vcf to disk
        uint32_t dir = pangraph_node_index / nbOfVCFsPerDir + 1; //get the good dir for this sample vcf
        auto vcfPath = VCFsDirs + "/" + int_to_string(dir) + "/" + prg_ptr->name + ".vcf";
        vcf.save(vcfPath, true, true, true, true, true, true, true);

        //add the vcf path to VCFPathsToBeConcatenated to concatenate after
        #pragma omp critical(VCFPathsToBeConcatenated)
        {
            VCFPathsToBeConcatenated.push_back(vcfPath);
        }


        if (genotype) {
            //genotype
            vcf.genotype(exp_depth_covgs, genotyping_error_rate, confidence_threshold, min_allele_covg_gt,
                         min_allele_fraction_covg_gt,
                         min_total_covg_gt, min_diff_covg_gt, false);

            //save the genotyped vcf to disk
            auto vcfGenotypedPath = VCFsGenotypedDirs + "/" + int_to_string(dir) + "/" + prg_ptr->name + "_genotyped.vcf";
            vcf.save(vcfGenotypedPath, true, true, true, true, true, true, true);

            //add the genotyped vcf path to VCFGenotypedPathsToBeConcatenated to concatenate after
            #pragma omp critical(VCFGenotypedPathsToBeConcatenated)
            {
                VCFGenotypedPathsToBeConcatenated.push_back(vcfGenotypedPath);
            }
        }
    }

    //generate all the multisample files
    vcf_ref_fa.save(outdir + "/pandora_multisample.vcf_ref.fa");
    VCF::concatenateVCFs(VCFPathsToBeConcatenated, outdir + "/pandora_multisample_consensus.vcf");
    if (genotype)
        VCF::concatenateVCFs(VCFGenotypedPathsToBeConcatenated, outdir + "/pandora_multisample_genotyped.vcf");


    // output a matrix/vcf which has the presence/absence of each prg in each sample
    BOOST_LOG_TRIVIAL(info) << "Output matrix";
    pangraph->save_matrix(outdir + "/pandora_multisample.matrix", sample_names);

    if (pangraph->nodes.empty()) {
        std::cout << "No LocalPRGs found to compare samples on. "
                  << "Is your genome_size accurate? Genome size is assumed to be "
                  << genome_size << " and can be updated with --genome_size" << std::endl;
    }

    // clear up
    index->clear();

    // current date/time based on current system
    std::cout << "FINISH: " << now() << std::endl;
    return 0;
}
