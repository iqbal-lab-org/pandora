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
using namespace std;

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

map<string, string> load_read_index(const string &readindex) {
    map<string, string> samples;
    string name, reads_path, line;
    ifstream myfile(readindex);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            istringstream linestream(line);
            if (std::getline(linestream, name, '\t')) {
                linestream >> reads_path;
                if (samples.find(name) != samples.end()) {
                    cout << "Warning: non-unique sample ids given! Only the last of these will be kept" << endl;
                }
                samples[name] = reads_path;
            }
        }
    } else {
        cerr << "Unable to open read index file " << readindex << endl;
        exit(1);
    }
    cout << now() << "Finished loading " << samples.size() << " samples from read index" << endl;
    return samples;
}

int pandora_compare(int argc, char *argv[]) {
    // if not enough arguments, print usage
    if (argc < 7) {
        show_compare_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile, readindex, outdir=".", vcf_refs_file;
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000, max_covg = 300; // default parameters
    int max_diff = 250;
    float e_rate = 0.11;
    bool illumina = false, clean = false, bin = false, genotype = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
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
        } else if((arg == "--max_covg")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_covg = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_covg option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotype")) {
            genotype = true;
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }

    //then run the programme...
    cout << "START: " << now() << endl;
    cout << "\nUsing parameters: " << endl;
    cout << "\tprgfile\t\t" << prgfile << endl;
    cout << "\treadindex\t" << readindex << endl;
    cout << "\toutdir\t" << outdir << endl;
    cout << "\tw\t\t" << w << endl;
    cout << "\tk\t\t" << k << endl;
    cout << "\tmax_diff\t" << max_diff << endl;
    cout << "\terror_rate\t" << e_rate << endl;
    cout << "\tvcf_refs\t" << vcf_refs_file << endl;
    cout << "\tillumina\t" << illumina << endl;
    cout << "\tclean\t" << clean << endl;
    cout << "\tbin\t" << bin << endl << endl;

    cout << now() << "Loading Index and LocalPRGs from file" << endl;
    Index *idx;
    idx = new Index();
    idx->load(prgfile, w, k);
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    // load read index
    cout << now() << "Loading read index file " << readindex << endl;
    map<string, string> samples = load_read_index(readindex);

    pangenome::Graph *pangraph, *pangraph_sample;
    pangraph = new pangenome::Graph();
    pangraph_sample = new pangenome::Graph();

    MinimizerHits *mhs;
    mhs = new MinimizerHits(100000);
    uint32_t covg;

    Fastaq consensus_fq(true, true);
    VCF master_vcf;

    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;
    vector<vector<uint32_t>> read_overlap_coordinates;

    // load vcf refs
    VCFRefs vcf_refs;
    string vcf_ref;
    if (!vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(vcf_refs_file, vcf_refs);
    }

    // for each sample, run pandora to get the sample pangraph
    for (map<string, string>::const_iterator sample = samples.begin(); sample != samples.end(); ++sample) {
        pangraph_sample->clear();
        mhs->clear();

        // make output dir for this sample
        string sample_outdir = outdir + "/" + sample->first;
        make_dir(sample_outdir + "/kmer_prgs");

        // construct the pangraph for this sample
        cout << now() << "Constructing pangenome::Graph from read file " << sample->second
             << " (this will take a while)" << endl;
        covg = pangraph_from_read_file(sample->second, mhs, pangraph_sample, idx, prgs, w, k, max_diff, e_rate,
                                       min_cluster_size, genome_size, illumina, clean, max_covg);

        cout << now() << "Finished with minihits, so clear " << endl;
        mhs->clear();

        cout << now() << "Writing pangenome::Graph to file " << sample_outdir << "/pandora.pangraph.gfa" << endl;
        write_pangraph_gfa(sample_outdir + "/pandora.pangraph.gfa", pangraph_sample);

        cout << now() << "Update LocalPRGs with hits" << endl;
        update_localPRGs_with_hits(pangraph_sample, prgs);

        cout << now() << "Estimate parameters for kmer graph model" << endl;
        estimate_parameters(pangraph_sample, sample_outdir, k, e_rate, covg, bin);

        cout << now() << "Find max likelihood PRG paths" << endl;
        for (auto c = pangraph_sample->nodes.begin(); c != pangraph_sample->nodes.end();) {
            if (!vcf_refs_file.empty()) {
                assert(vcf_refs.find(prgs[c->second->prg_id]->name) != vcf_refs.end());
                vcf_ref = vcf_refs[prgs[c->second->prg_id]->name];
            }

            prgs[c->second->prg_id]->add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp, w, bin, covg);
            if (kmp.empty())
            {
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
    }

    // for each pannode in graph, find a best reference and output a vcf and aligned fasta of sample paths through it
    cout << now() << "Multi-sample pangraph has " << pangraph->nodes.size() << " nodes" << endl;
    for (const auto &c: pangraph->nodes) {
        cout << " c.first: " << c.first;
        cout << " prgs[c.first]->name: " << prgs[c.first]->name << endl;

        if (!vcf_refs_file.empty()
            and vcf_refs.find(prgs[c.second->prg_id]->name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prgs[c.second->prg_id]->name];
        }

        string node_outdir = outdir + "/" + c.second->get_name();
        make_dir(node_outdir);

        c.second->output_samples(prgs[c.first], node_outdir, w, vcf_ref);
    }
    if(genotype) {
        master_vcf.genotype(covg,0.01,30,false);
        master_vcf.save(outdir + "/pandora_genotyped.vcf" , true, true, true, true, false, false, false);
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    cout << now() << "Output matrix" << endl;
    pangraph->save_matrix(outdir + "multisample.matrix");

    // clear up
    cout << now() << "Clear up" << endl;
    for (uint32_t j = 0; j != prgs.size(); ++j) {
        delete prgs[j];
    }
    idx->clear();
    delete idx;
    mhs->clear();
    delete mhs;
    pangraph->clear();
    delete pangraph;
    pangraph_sample->clear();
    delete pangraph_sample;

    // current date/time based on current system
    cout << "FINISH: " << now() << endl;
    return 0;
}
