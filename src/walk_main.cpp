#include <cstring>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include "localPRG.h"
#include "utils.h"
#include "localgraph.h"
#include "localnode.h"
#include "fastaq_handler.h"

using namespace std;

int pandora_walk(int argc, char *argv[]) // the "pandora walk" comand
{
    if (argc != 3) {
        fprintf(stderr, "Usage: pandora walk <in_prg.fa> [<seq.fa> | --top | --bottom]\n");
        return 1;
    }

    // load prgs
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, argv[1]);

    vector<LocalNodePtr> npath;

    if (strcmp(argv[2], "--top") == 0) {
        for (auto prg_ptr : prgs){
            npath = prg_ptr->prg.top_path();
            cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                cout << "->" << npath[j]->id;
            }
            cout << endl;
        }
        return 0;
    } else if (strcmp(argv[2], "--bottom") == 0) {
        for (auto prg_ptr : prgs){
            npath = prg_ptr->prg.bottom_path();
            cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                cout << "->" << npath[j]->id;
            }
            cout << endl;
        }
        return 0;
    }

    // for each read in readfile,  infer node path along sequence
    string read_string;
    FastaqHandler readfile(argv[2]);
    while(not readfile.eof()) {
        readfile.get_next();
        //cout << "Try to find gene " << readfile.num_reads_parsed << " " << readfile.name << endl << readfile.read << endl;
        for (auto prg_ptr : prgs){
            npath = prg_ptr->prg.nodes_along_string(readfile.read);
            if (not npath.empty()) {
                cout << readfile.name << "\t" << prg_ptr->name << "\t";
                for (uint32_t j = 0; j != npath.size(); ++j) {
                    cout << "->" << npath[j]->id;
                }
                cout << endl;
            }
        }
    }

    return 0;
}

