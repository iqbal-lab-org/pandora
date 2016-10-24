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
#include "pangraph.h"
#include "index.h"

using std::set;
using std::vector;
using namespace std;

int main(int argc, char *argv[])
{
    assert(argc >= 3 && "Not enough input arguments - should have PRG file, read file, and an output prefix!");
    // current date/time based on current system
    time_t now = time(0);
    // convert now to string form
    char* dt = ctime(&now);
    cout << "START: " << dt << endl;

    uint32_t w = 1;//5;
    uint32_t k = 3;//15;
    int max_diff = 1;//500;
    uint32_t cluster_thresh = 1;//4;

    cout << "Building Index from PRG file" << endl;
    Index *idx;
    idx = new Index();
    vector<LocalPRG*> prgs;
    index_prg_file(prgs, argv[1], idx, w, k);

    cout << "Reading read file" << endl;
    PanGraph *pangraph;
    pangraph = new PanGraph();
    pangraph_from_read_file(argv[2], pangraph, idx, w, k, max_diff, cluster_thresh);
    cout << "Finished building pangraph from reads" << endl;
	
    // delete all the Seq objects
    for (uint32_t j=0; j<prgs.size(); ++j)
    {
	delete prgs[j];
    }

    // current date/time based on current system
    now = time(0);
    dt = ctime(&now);
    cout << "FINISH: " << dt << endl;
    return 0;
}

