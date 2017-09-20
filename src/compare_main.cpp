/*
 * C++ Program to call variation between a set of samples. 
 * To do this, runs pandora map process on each sample/readfile to find the best gene sequence for 
 * each gene inferred present in that sample. Then compares the paths found through each gene and 
 * outputs a vcf and aligned fasta for each one. 
 */
/*
 * QUESTIONS:
 * How do I handle multiple occurances of a gene in a sample? I would like to output all the versions called against my new reference
 */

#include <iostream>
#include <sstream>
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
#include "pangraph.h"
#include "pannode.h"
#include "index.h"
#include "estimate_parameters.h"

using std::set;
using std::vector;
using namespace std;

static void show_compare_usage()
{
    std::cerr << "Usage: pandora compare -p PRG_FILE -r READ_INDEX -o OUT_PREFIX <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-p,--prg_file PRG_FILE\t\tSpecify a fasta-style prg file\n"
	      << "\t-r,--read_index READ_INDEX\tSpecify a file with a line per sample\n"
	      << "\t\t\t\t\tsample_id <tab> filepath to reads in fasta/q format\n"
	      << "\t-o,--out_prefix OUT_PREFIX\tSpecify prefix of output\n"
	      << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, default 1\n"
	      << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
	      << "\t-m,--max_diff INT\t\tMaximum distance between consecutive hits within a cluster, default 500 (bps)\n"
	      << "\t-e,--error_rate FLOAT\t\tEstimated error rate for reads, default 0.11\n"
	      << "\t--output_kg\t\t\tSave kmer graphs with fwd and rev coverage annotations for found localPRGs\n"
	      << "\t--output_vcf\t\t\tSave a vcf file for each found localPRG\n"
	      << "\t--method\t\t\tMethod for path inference, can be max likelihood (default), 'min' to maximize\n"
	      << "\t\t\t\t\tthe min probability on the path, or 'both' to create outputs with both methods\n"
	      << "\t--output_comparison_paths\tSave a fasta file for a random selection of paths through localPRG\n"
              << std::endl;
}

map<string,string> load_read_index(const string& readindex)
{
    map<string,string> samples;
    string name, reads_path, line;
    ifstream myfile (readindex);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            istringstream linestream(line);
            if (std::getline(linestream, name, '\t'))
            {
                linestream >> reads_path;
                if (samples.find(name) != samples.end())
                {
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

int pandora_compare(int argc, char* argv[])
{
    // if not enough arguments, print usage
    if (argc < 7) {
        show_compare_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile, readindex, prefix, readfile;
    uint32_t w=1, k=15; // default parameters
    int max_diff = 500;
    float e_rate = 0.11;
    bool output_kg=false, output_vcf=false, max_path=true, min_path=false, output_comparison_paths=false;
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
                readfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--read_index option requires one argument." << std::endl;
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
                w = (unsigned)atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "-w option requires one argument." << std::endl;
                return 1;
            }  
	} else if (arg == "-k") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                k = (unsigned)atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
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
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--error_rate option requires one argument." << std::endl;
                return 1;
            }
	} else if ((arg == "--output_kg")) {
	    output_kg = true;
	} else if ((arg == "--output_vcf")) {
            output_vcf = true;
	} else if ((arg == "--method")) {
	    if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                string method = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
		if (method == "min")
		{
		    max_path = false;
		    min_path = true;
		} else if (method == "both")
		{
		    min_path = true;
		}
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--method option requires one argument." << std::endl;
                return 1;
            }
	} else if ((arg == "--output_comparison_paths")) {
            output_comparison_paths = true;
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }

    //then run the programme...
    cout << "START: " << now() << endl;
    cout << "\nUsing parameters: " << endl;
    cout << "\tprgfile\t\t" << prgfile << endl;
    cout << "\treadindex\t" << readindex << endl;
    cout << "\tout_prefix\t" << prefix << endl;
    cout << "\tw\t\t" << w << endl;
    cout << "\tk\t\t" << k << endl;
    cout << "\tmax_diff\t" << max_diff << endl;
    cout << "\terror_rate\t" << e_rate << endl;
    cout << "\toutput_kg\t" << output_kg << endl;
    cout << "\toutput_vcf\t" << output_vcf << endl;
    cout << "\tmax_path\t" << max_path << endl;
    cout << "\tmin_path\t" << min_path << endl;
    cout << "\toutput_comparison_paths\t" << output_comparison_paths << endl << endl;

    cout << now() << "Loading Index and LocalPRGs from file" << endl;
    Index *idx;
    idx = new Index();
    idx->load(prgfile, w, k);
    vector<LocalPRG*> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    PanGraph *pangraph, *pangraph_sample;
    pangraph = new PanGraph();
    pangraph_sample = new PanGraph();
    MinimizerHits *mhs;
    mhs = new MinimizerHits(100*idx->minhash.size());

    // load read index
    map<string,string> samples = load_read_index(readindex);
    vector<KmerNode*> kmp;

    // for each sample, run pandora to get the sample pangraph
    for (map<string,string>::const_iterator sample = samples.begin(); sample!=samples.end(); ++sample)
    {
	pangraph_sample->clear();
	mhs->clear();
	
	// construct the pangraph for this sample
        cout << now() << "Constructing PanGraph from read file " << sample->second << endl;
        pangraph_from_read_file(sample->second, mhs, pangraph_sample, idx, prgs, w, k, max_diff);
    
        cout << now() << "Update LocalPRGs with hits" << endl;
        update_localPRGs_with_hits(pangraph_sample, prgs);

        cout << now() << "Estimate parameters for kmer graph model" << endl;
        estimate_parameters(pangraph_sample, prefix, k, e_rate);

        cout << now() << "Find max likelihood PRG paths and write to files:" << endl;
        for (auto c: pangraph_sample->nodes)
        {
	    kmp = prgs[c.second->prg_id]->find_path_and_variants(c.second, prefix, w, max_path, min_path, output_vcf, output_comparison_paths);
	    pangraph->add_node(c.second->prg_id, c.second->name, sample->first, kmp, prgs[c.second->prg_id]);
        }
    }

    // for each pannode in graph, find a best reference and output a vcf and aligned fasta of sample paths through it
    for (auto c: pangraph->nodes)
    {
	c.second->output_samples_vcf(prgs[c.first], prefix, w);
    }

    // output a matrix/vcf which has the presence/absence of each prg in each sample
    pangraph->save_matrix(prefix + "multisample.matrix");

    // clear up
    for (uint32_t j=0; j!=prgs.size(); ++j)
    {
	delete prgs[j];
    }
    idx->clear();
    delete idx;
    delete mhs;
    delete pangraph;
    delete pangraph_sample;

    // current date/time based on current system
    cout << "FINISH: " << now() << endl;
    return 0;
}
