/*
 * C++ Program to find (w,k)-minimizer hits between a fasta of reads and a set of gene PRGs and output a gfa of genes that there is evidence have been traversed and the order they have been seen.
 */
/*
 * To do:
 * 1. Allow reverse complement strand minimizers
 * 2. Make sure mapped region not got large indels from unsupported intervals
 * 3. Command line argument handling for option of dumped_db screening
 * 4. Change structure so not copying minimzers etc, but use pointers
 */
#include <iostream>
#include <ctime>
#include <cstring>
#include <vector>
#include <set>
#include <tuple>
#include <functional>
#include <ctype.h>
#include <fstream>
#include <algorithm>
#include <map>
#include <assert.h>

#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "estimate_parameters.h"
#include "noise_filtering.h"

using std::set;
using std::vector;
using namespace std;

static void show_map_usage() {
    std::cerr << "Usage: pandora map -p PRG_FILE -r READ_FILE -o OUT_PREFIX <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-p,--prg_file PRG_FILE\t\tSpecify a fasta-style prg file\n"
              << "\t-r,--read_file READ_FILE\tSpecify a file of reads in fasta format\n"
              << "\t-o,--out_prefix OUT_PREFIX\tSpecify prefix of output\n"
              << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, must be <=k, default 14\n"
              << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
              << "\t-m,--max_diff INT\t\tMaximum distance between consecutive hits within a cluster, default 500 (bps)\n"
              << "\t-e,--error_rate FLOAT\t\tEstimated error rate for reads, default 0.11\n"
              << "\t--genome_size\tNUM_BP\tEstimated length of genome, used for coverage estimation\n"
              << "\t--output_kg\t\t\tSave kmer graphs with fwd and rev coverage annotations for found localPRGs\n"
              << "\t--output_vcf\t\t\tSave a vcf file for each found localPRG\n"
              << "\t--vcf_refs REF_FASTA\t\tA fasta file with an entry for each LocalPRG giving reference sequence for\n"
              << "\t\t\t\t\tVCF. Must have a perfect match in the graph and the same name as the graph\n"
              << "\t--output_comparison_paths\tSave a fasta file for a random selection of paths through localPRG\n"
              << "\t--output_covgs\tSave a file of covgs for each localPRG present, one number per base of fasta file\n"
              << "\t--illumina\t\t\tData is from illumina rather than nanopore, so is shorter with low error rate\n"
              << "\t--clean\t\t\tAdd a step to clean and detangle the pangraph\n"
              << "\t--nbin\t\t\tUse negative binomial model for kmer coverages\n"
              << std::endl;
}

int pandora_map(int argc, char *argv[]) {
    // if not enough arguments, print usage
    if (argc < 7) {
        show_map_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile, readfile, prefix, vcf_refs_file;
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000; // default parameters
    int max_diff = 250;
    float e_rate = 0.11;
    bool output_kg = false, output_vcf = false;
    bool output_comparison_paths = false, illumina = false, clean = false;
    bool output_covgs = false, nbin = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_map_usage();
            return 0;
        } else if ((arg == "-p") || (arg == "--prg_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prgfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--prg_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-r") || (arg == "--read_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                readfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--read_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-o") || (arg == "--out_prefix")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prefix = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--out_prefix option requires one argument." << std::endl;
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
                e_rate = static_cast<float>(atof(
                        argv[++i])); // Increment 'i' so we don't get the argument as the next argv[i].
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
        } else if ((arg == "--output_kg")) {
            output_kg = true;
        } else if ((arg == "--output_vcf")) {
            output_vcf = true;
        } else if (arg == "--vcf_refs") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                vcf_refs_file = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--vcf_refs option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--output_covgs")) {
            output_covgs = true;
        } else if ((arg == "--output_comparison_paths")) {
            output_comparison_paths = true;
        } else if ((arg == "--illumina")) {
            illumina = true;
            if (e_rate > 0.05) {
                e_rate = 0.001;
            }
        } else if ((arg == "--clean")) {
            clean = true;
        } else if ((arg == "--nbin")) {
            nbin = true;
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }

    assert(w <= k);
    assert(not prgfile.empty());

    //then run the programme...
    cout << "START: " << now() << endl;
    cout << "\nUsing parameters: " << endl;
    cout << "\tprgfile\t\t" << prgfile << endl;
    cout << "\treadfile\t" << readfile << endl;
    cout << "\tout_prefix\t" << prefix << endl;
    cout << "\tw\t\t" << w << endl;
    cout << "\tk\t\t" << k << endl;
    cout << "\tmax_diff\t" << max_diff << endl;
    cout << "\terror_rate\t" << e_rate << endl;
    cout << "\toutput_kg\t" << output_kg << endl;
    cout << "\toutput_vcf\t" << output_vcf << endl;
    cout << "\tvcf_refs\t" << vcf_refs_file << endl;
    cout << "\toutput_comparison_paths\t" << output_comparison_paths << endl;
    cout << "\toutput_covgs\t" << output_covgs << endl;
    cout << "\tillumina\t" << illumina << endl;
    cout << "\tclean\t" << clean << endl;
    cout << "\tnbin\t" << nbin << endl << endl;

    cout << now() << "Loading Index and LocalPRGs from file" << endl;
    Index *idx;
    idx = new Index();
    idx->load(prgfile, w, k);
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    cout << now() << "Constructing pangenome::Graph from read file (this will take a while)" << endl;
    MinimizerHits *mhs;
    mhs = new MinimizerHits(100 * idx->minhash.size());
    pangenome::Graph *pangraph;
    pangraph = new pangenome::Graph();
    uint covg = pangraph_from_read_file(readfile, mhs, pangraph, idx, prgs, w, k, max_diff, e_rate, min_cluster_size,
                                        genome_size, illumina, clean);

    cout << now() << "Writing pangenome::Graph to file " << prefix << ".pangraph.gfa" << endl;
    write_pangraph_gfa(prefix + ".pangraph.gfa", pangraph);

    cout << now() << "Update LocalPRGs with hits" << endl;
    update_localPRGs_with_hits(pangraph, prgs);

    cout << now() << "Estimate parameters for kmer graph model" << endl;
    estimate_parameters(pangraph, prefix, k, e_rate, covg, nbin);

    cout << now() << "Find PRG paths and write to files:" << endl;
    VCFRefs vcf_refs;
    string vcf_ref;
    vector<KmerNodePtr> kmp;
    if (output_vcf and !vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(vcf_refs_file, vcf_refs);
    }
    for (auto c: pangraph->nodes) {
        if (output_vcf
            and !vcf_refs_file.empty()
            and vcf_refs.find(prgs[c.second->prg_id]->name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prgs[c.second->prg_id]->name];
        }
        kmp = prgs[c.second->prg_id]->find_path_and_variants(c.second, prefix, w, vcf_ref, output_vcf,
                                                       output_comparison_paths, output_covgs, nbin);
        if (kmp.empty())
        {
            continue;
            //could also remove the node from pg
        }

        if (output_kg) {
            c.second->kmer_prg.save(prefix + "." + c.second->get_name() + ".kg.gfa", prgs[c.second->prg_id]);
        }
        //prgs[c.second->id]->kmer_prg.save_covg_dist(prefix + "." + prgs[c.second->id]->name + ".covg.txt");
        //cout << "\t\t" << prefix << "." << prgs[c.second->id]->name << ".gfa" << endl;
        //prgs[c.second->id]->prg.write_gfa(prefix + "." + prgs[c.second->id]->name + ".gfa");
    }

    //cout << now() << "Writing LocalGraphs to files:" << endl;	
    // for each found localPRG, also write out a gfa 
    // then delete the localPRG object
    for (uint32_t j = 0; j != prgs.size(); ++j) {
        //cout << "\t\t" << prefix << "_" << prgs[j]->name << ".gfa" << endl;
        //prgs[j]->prg.write_gfa(prefix + "_" + prgs[j]->name + ".gfa");
        delete prgs[j];
    }
    idx->clear();
    delete idx;
    mhs->clear();
    delete mhs;
    pangraph->clear();
    delete pangraph;

    // current date/time based on current system
    cout << "FINISH: " << now() << endl;
    return 0;
}

