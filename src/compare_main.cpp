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
              << "\t--genome_size\tNUM_BP\tEstimated length of genome, used for coverage estimation\n"
              << "\t--vcf_refs REF_FASTA\t\tA fasta file with an entry for each LocalPRG giving reference sequence for\n"
              << "\t\t\t\t\tVCF. Must have a perfect match in the graph and the same name as the graph\n"
              << "\t--illumina\t\t\tData is from illumina rather than nanopore, so is shorter with low error rate\n"
              << "\t--clean\t\t\tAdd a step to clean and detangle the pangraph\n"
              << "\t--bin\t\t\tUse binomial model for kmer coverages, default is negative binomial\n"
              << "\t--max_covg\t\t\tMaximum average coverage from reads to accept\n"
              << "\t--genotype\t\t\tAdd extra step to carefully genotype sites\n"
              << std::endl;
}

using SampleIdText = std::string;
using SampleFpath = std::string;

std::map<SampleIdText, SampleFpath> load_read_index(const std::string &read_index_fpath) {
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
    return samples;
}

int pandora_compare(int argc, char *argv[]) {
    // if not enough arguments, print usage
    if (argc < 7) {
        show_compare_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    std::string prgfile, read_index_fpath, outdir = ".", vcf_refs_file;
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000, max_covg = 300; // default parameters
    int max_diff = 250;
    float e_rate = 0.11;
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
                w = (unsigned) atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-w option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-k") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                k = (unsigned) atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-k option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-m") || (arg == "--max_diff")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_diff = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_diff option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-e") || (arg == "--error_rate")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                e_rate = atof(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
                if (e_rate < 0.01) {
                    illumina = true;
                }
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--error_rate option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genome_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                genome_size = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
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
            if (e_rate > 0.05) {
                e_rate = 0.001;
            }
        } else if ((arg == "--clean")) {
            clean = true;
        } else if ((arg == "--bin")) {
            bin = true;
        } else if ((arg == "--max_covg")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_covg = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_covg option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotype")) {
            genotype = true;
        } else {
            std::cerr << argv[i] << " could not be attributed to any parameter" << std::endl;
        }
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
    std::cout << "\tvcf_refs\t" << vcf_refs_file << std::endl;
    std::cout << "\tillumina\t" << illumina << std::endl;
    std::cout << "\tclean\t" << clean << std::endl;
    std::cout << "\tbin\t" << bin << std::endl << std::endl;

    std::cout << now() << "Loading Index and LocalPRGs from file" << std::endl;
    auto index = std::make_shared<Index>();
    index->load(prgfile, w, k);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    // load read index
    std::cout << now() << "Loading read index file " << read_index_fpath << std::endl;
    auto samples = load_read_index(read_index_fpath);

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    auto pangraph_sample = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits(100000));

    Fastaq consensus_fq(true, true);

    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;
    vector<vector<uint32_t>> read_overlap_coordinates;

    // for each sample, run pandora to get the sample pangraph
    uint32_t last_covg = 0;
    uint32_t sample_id = 0;
    for (const auto &sample: samples) {
        pangraph_sample->clear();
        minimizer_hits->clear();

        const auto &sample_name = sample.first;
        const auto &sample_fpath = sample.second;

        // make output dir for this sample
        auto sample_outdir = outdir + "/" + sample_name;
        fs::create_directories(sample_outdir);

        // construct the pangraph for this sample
        BOOST_LOG_TRIVIAL(info) << "Constructing pangenome::Graph from read file "
                                << sample_fpath
                                << " (this will take a while)";
        uint32_t covg = pangraph_from_read_file(sample_fpath,
                                                minimizer_hits,
                                                pangraph_sample,
                                                index, prgs, w, k,
                                                max_diff, e_rate,
                                                min_cluster_size, genome_size, illumina, clean, max_covg);
        last_covg = covg;
        BOOST_LOG_TRIVIAL(info) << "Finished with minihits, so clear ";
        minimizer_hits->clear();

        BOOST_LOG_TRIVIAL(info) << "Writing pangenome::Graph to file "
                                << sample_outdir << "/pandora.pangraph.gfa";
        write_pangraph_gfa(sample_outdir + "/pandora.pangraph.gfa", pangraph_sample);

        if (pangraph_sample->nodes.empty()) {
            BOOST_LOG_TRIVIAL(warning) << "Found no LocalPRGs in the reads for sample " << sample_name;
        }

        BOOST_LOG_TRIVIAL(info) << "Update LocalPRGs with hits";
        pangraph_sample->setup_kmergraphs(prgs, 1);
        pangraph_sample->add_hits_to_kmergraphs(prgs, 0);

        BOOST_LOG_TRIVIAL(info) << "Estimate parameters for kmer graph model";
        estimate_parameters(pangraph_sample, sample_outdir, k, e_rate, covg, bin, 0);

        BOOST_LOG_TRIVIAL(info) << "Find max likelihood PRG paths";
        auto sample_pangraph_size = pangraph_sample->nodes.size();
        for (auto c = pangraph_sample->nodes.begin(); c != pangraph_sample->nodes.end();) {
            LocalPRG local_prg = *prgs[c->second->prg_id];
            local_prg.add_consensus_path_to_fastaq(consensus_fq,
                                                   c->second,
                                                   kmp, lmp, w,
                                                   bin, covg, 0);

            if (kmp.empty()) {
                c = pangraph_sample->remove_node(c->second);
                continue;
            }
            pangraph->add_node(c->second->prg_id, c->second->name, sample_name, sample_id,
                               prgs[c->second->prg_id], kmp);

            ++c;
        }

        pangraph->setup_kmergraphs(prgs, samples.size());
        pangraph->copy_coverages_to_kmergraphs(*pangraph_sample, sample_id);

        consensus_fq.save(sample_outdir + "/pandora.consensus.fq.gz");
        consensus_fq.clear();
        if (pangraph_sample->nodes.empty() and sample_pangraph_size > 0) {
            std::cout << "WARNING: All LocalPRGs found were removed for sample " << sample_name
                      << ". Is your genome_size accurate? Genome size is assumed to be " << genome_size
                      << " and can be updated with --genome_size" << std::endl;
        }

        sample_id++;
    }

    // for each pannode in graph, find a best reference
    // and output a vcf and aligned fasta of sample paths through it
    BOOST_LOG_TRIVIAL(info) << "Multi-sample pangraph has "
                            << pangraph->nodes.size() << " nodes";

    // load vcf refs
    VCFRefs vcf_refs;
    std::string vcf_ref;
    if (!vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(vcf_refs_file, vcf_refs);
    }

    VCF master_vcf;

    for (const auto &pangraph_node_entry: pangraph->nodes) {
        BOOST_LOG_TRIVIAL(debug) << "Consider next node";
        const auto &node_id = pangraph_node_entry.first;
        pangenome::Node &pangraph_node = *pangraph_node_entry.second;

        const auto &prg_id = pangraph_node.prg_id;
        assert(prgs.size() > prg_id);
        const auto& prg_ptr = prgs[prg_id];

        const auto vcf_reference_path = pangraph->infer_node_vcf_reference_path(pangraph_node, prg_ptr, w, vcf_refs);
        BOOST_LOG_TRIVIAL(debug) << " c.first: " << node_id << " prgs[c.first]->name: " << prg_ptr->name;

        pangraph_node.construct_multisample_vcf(master_vcf, vcf_reference_path, prg_ptr, w);
    }
    master_vcf.save(outdir + "/pandora_multisample_consensus.vcf", true, true, true, true, true, true, true);

    if (genotype) {
        master_vcf.genotype(last_covg, 0.01, 30, false);
        master_vcf.save(outdir + "/pandora_multisample_genotyped.vcf", true, true, true, true, true, true, true);
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    BOOST_LOG_TRIVIAL(info) << "Output matrix";
    pangraph->save_matrix(outdir + "/pandora_multisample.matrix");

    if (pangraph->nodes.empty()) {
        std::cout << "No LocalPRGs found to compare samples on. "
                  << "Is your genome_size accurate? Genome size is assumed to be "
                  << genome_size << " and can be updated with --genome_size" << std::endl;
    }

    // clear up
    index->clear();
    minimizer_hits->clear();
    pangraph->clear();
    pangraph_sample->clear();

    // current date/time based on current system
    std::cout << "FINISH: " << now() << std::endl;
    return 0;
}
