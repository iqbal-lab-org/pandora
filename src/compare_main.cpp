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

std::map<std::string, std::string> load_read_index(const std::string &readindex) {
    std::map<std::string, std::string> samples;
    std::string name, reads_path, line;
    std::ifstream myfile(readindex);
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
        std::cerr << "Unable to open read index file " << readindex << std::endl;
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
    std::string prgfile, readindex, outdir = ".", vcf_refs_file;
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
                readindex = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
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
    std::cout << "\treadindex\t" << readindex << std::endl;
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
    std::cout << now() << "Loading read index file " << readindex << std::endl;
    auto samples = load_read_index(readindex);

    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    auto pangraph_sample = std::make_shared<pangenome::Graph>(pangenome::Graph());

    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits(100000));
    uint32_t covg;

    Fastaq consensus_fq(true, true);
    VCF master_vcf;

    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;
    vector<vector<uint32_t>> read_overlap_coordinates;

    // load vcf refs
    VCFRefs vcf_refs;
    std::string vcf_ref;
    if (!vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(vcf_refs_file, vcf_refs);
    }

    // for each sample, run pandora to get the sample pangraph
    for (auto sample = samples.begin(); sample != samples.end(); ++sample) {
        pangraph_sample->clear();
        minimizer_hits->clear();

        // make output dir for this sample
        auto sample_outdir = outdir + "/" + sample->first;
        fs::create_directories(sample_outdir + "/kmer_prgs");

        // construct the pangraph for this sample
        std::cout << now() << "Constructing pangenome::Graph from read file " << sample->second
                  << " (this will take a while)" << std::endl;
        covg = pangraph_from_read_file(sample->second, minimizer_hits, pangraph_sample, index, prgs, w, k, max_diff, e_rate,
                                       min_cluster_size, genome_size, illumina, clean, max_covg);

        std::cout << now() << "Finished with minihits, so clear " << std::endl;
        minimizer_hits->clear();

        std::cout << now() << "Writing pangenome::Graph to file " << sample_outdir << "/pandora.pangraph.gfa"
                  << std::endl;
        write_pangraph_gfa(sample_outdir + "/pandora.pangraph.gfa", pangraph_sample);

        if (pangraph_sample->nodes.empty()) {
            std::cout << "WARNING: Found no LocalPRGs in the reads for sample " << sample->first << std::endl;
        }

        std::cout << now() << "Update LocalPRGs with hits" << std::endl;
        update_localPRGs_with_hits(pangraph_sample, prgs);

        std::cout << now() << "Estimate parameters for kmer graph model" << std::endl;
        estimate_parameters(pangraph_sample, sample_outdir, k, e_rate, covg, bin);

        std::cout << now() << "Find max likelihood PRG paths" << std::endl;
        auto sample_pangraph_size = pangraph_sample->nodes.size();
        for (auto c = pangraph_sample->nodes.begin(); c != pangraph_sample->nodes.end();) {
            if (!vcf_refs_file.empty()) {
                assert(vcf_refs.find(prgs[c->second->prg_id]->name) != vcf_refs.end());
                vcf_ref = vcf_refs[prgs[c->second->prg_id]->name];
            }

            prgs[c->second->prg_id]->add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp, w, bin, covg);
            if (kmp.empty()) {
                c = pangraph_sample->remove_node(c->second);
                continue;
            }

            prgs[c->second->prg_id]->add_variants_to_vcf(master_vcf, c->second, vcf_ref, kmp, lmp);

            pangraph_sample->save_kmergraph_coverages(outdir, c->second->name);
            pangraph->add_node(c->second->prg_id, c->second->name, sample->first, kmp, prgs[c->second->prg_id]);

            ++c;
        }
        consensus_fq.save(sample_outdir + "/pandora.consensus.fq.gz");
        consensus_fq.clear();
        //master_vcf.save(sample_outdir + "/pandora_consensus.vcf" , true, true, true, true, true, true, true);
        //master_vcf.clear();
        if (pangraph_sample->nodes.empty() and sample_pangraph_size > 0) {
            std::cout << "WARNING: All LocalPRGs found were removed for sample " << sample->first
                      << ". Is your genome_size accurate? Genome size is assumed to be " << genome_size
                      << " and can be updated with --genome_size" << std::endl;
        }
    }

    // for each pannode in graph, find a best reference and output a vcf and aligned fasta of sample paths through it
    std::cout << now() << "Multi-sample pangraph has " << pangraph->nodes.size() << " nodes" << std::endl;
    for (const auto &c: pangraph->nodes) {
        std::cout << " c.first: " << c.first;
        std::cout << " prgs[c.first]->name: " << prgs[c.first]->name << std::endl;

        if (!vcf_refs_file.empty()
            and vcf_refs.find(prgs[c.second->prg_id]->name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prgs[c.second->prg_id]->name];
        }

        auto node_outdir = outdir + "/" + c.second->get_name();
        fs::create_directories(node_outdir);

        c.second->output_samples(prgs[c.first], node_outdir, w, vcf_ref);
    }
    if (genotype) {
        master_vcf.genotype(covg, 0.01, 30, false);
        master_vcf.save(outdir + "/pandora_genotyped.vcf", true, true, true, true, false, false, false);
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    std::cout << now() << "Output matrix" << std::endl;
    pangraph->save_matrix(outdir + "/multisample.matrix");

    // clear up
    std::cout << now() << "Clear up" << std::endl;
    index->clear();
    minimizer_hits->clear();
    pangraph->clear();
    pangraph_sample->clear();

    if (pangraph->nodes.empty()) {
        std::cout << "No LocalPRGs found to compare samples on. "
                  << "Is your genome_size accurate? Genome size is assumed to be "
                  << genome_size << " and can be updated with --genome_size" << std::endl;
    }

    // current date/time based on current system
    std::cout << "FINISH: " << now() << std::endl;
    return 0;
}
