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
#include "pangraph.h"
#include "index.h"

using std::set;
using std::vector;
using namespace std;

static void show_usage(string name)
{
    std::cerr << "Usage: " << name << " -p PRG_FILE -r READ_FILE -o OUT_PREFIX <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-p,--prg_file PRG_FILE\t\tSpecify a fasta-style prg file\n"
	      << "\t-r,--read_file READ_FILE\tSpecify a file of reads in fasta format\n"
	      << "\t-o,--out_prefix OUT_PREFIX\tSpecify prefix of output\n"
	      << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers\n"
	      << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers\n"
	      << "\t-m,--max_diff MAX_DIFF\t\tMaximum distance between consecutive hits within a cluster\n"
	      << "\t-c,--cthresh CLUSTER_THRESH\tMinimum number of hits in a cluster"
              << std::endl;
}

int main(int argc, char* argv[])
{
    // if not enough arguments, print usage
    if (argc < 7) {
        show_usage(argv[0]);
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile, readfile, prefix;
    uint32_t w=1, k=3, cluster_thresh = 1; // default parameters
    int max_diff = 1;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else if ((arg == "-p") || (arg == "--prg_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prgfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
		cout << "prgfile: " << prgfile << endl;
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--prg_file option requires one argument." << std::endl;
                return 1;
            }
	} else if ((arg == "-r") || (arg == "--read_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                readfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
	        cout << "readfile: " << readfile << endl;
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--read_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-o") || (arg == "--out_prefix")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prefix = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
		cout << "prefix: " << prefix << endl;
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
	} else if ((arg == "-c") || (arg == "--cthresh")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                cluster_thresh = (unsigned)atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                  std::cerr << "--cthresh option requires one argument." << std::endl;
                return 1;
            }
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }

    //then run the programme...
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    string dt = ctime(&now);
    cout << "START: " << dt << endl;

    now = time(0);
    dt = ctime(&now);
    string sdt = dt.substr(0,dt.length()-1);
    cout << sdt << " Building Index from PRG file" << endl;
    Index *idx;
    idx = new Index();
    vector<LocalPRG*> prgs;
    index_prg_file(prgs, prgfile, idx, w, k);

    now = time(0);
    dt = ctime(&now);
    sdt = dt.substr(0,dt.length()-1);
    cout << sdt << " Constructing PanGraph from read file" << endl;
    PanGraph *pangraph;
    pangraph = new PanGraph();
    pangraph_from_read_file(readfile, pangraph, idx, prgs, w, k, max_diff, cluster_thresh);

    now = time(0);
    dt = ctime(&now);
    sdt = dt.substr(0,dt.length()-1);
    cout << sdt << " Writing PanGraph to file " << prefix << "_pangraph.gfa" << endl;
    pangraph->write_gfa(prefix + "_pangraph.gfa");

    now = time(0);
    dt = ctime(&now);
    sdt = dt.substr(0,dt.length()-1);
    cout << sdt << " Writing LocalGraphs to files:" << endl;	
    // for each found localPRG, also write out a gfa 
    // then delete the localPRG object
    for (uint32_t j=0; j<prgs.size(); ++j)
    {
        cout << "\t\t" << prefix << "_" << prgs[j]->name << ".gfa" << endl;
	prgs[j]->prg.write_gfa(prefix + "_" + prgs[j]->name + ".gfa");
	delete prgs[j];
    }
    delete idx;
    delete pangraph;

    // current date/time based on current system
    now = time(0);
    dt = ctime(&now);
    cout << "FINISH: " << dt << endl;
    return 0;
}

