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

int pandora_get_vcf_ref(int argc, char *argv[]) // the "pandora walk" comand
{
    if (argc != 3) {
        fprintf(stderr, "Usage: pandora get_vcf_ref <in_prg.fa> <seq.fa>\n");
        return 1;
    }

    // load prgs
    vector<LocalPRG *> prgs;
    read_prg_file(prgs, argv[1]);

    vector<LocalNodePtr> npath;
    string read_string;
    FastaqHandler readfile(argv[2]);

    for (auto prg_ptr : prgs){
        readfile.get_id(0);
        while(not readfile.eof()) {
            npath = prg_ptr->prg.nodes_along_string(readfile.read);
            if (not npath.empty()) {
                cout << ">" << prg_ptr->name << endl << readfile.read << endl;
                break;
            }
            readfile.get_next();
        }
        if (not npath.empty()) {
            cout << ">" << prg_ptr->name << endl << readfile.read << endl;
        } else {
            npath = prg_ptr->prg.top_path();
            cout << ">" << prg_ptr->name << endl << prg_ptr->string_along_path(npath) << endl;
        }
    }

    return 0;
}

