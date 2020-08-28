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

int pandora_walk(int argc, char* argv[]) // the "pandora walk" comand
{
    if (argc != 3) {
        fprintf(
            stderr, "Usage: pandora walk <in_prg.fa> [<seq.fa> | --top | --bottom]\n");
        return 1;
    }

    // load prgs
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, argv[1]);

    std::vector<LocalNodePtr> npath;

    if (strcmp(argv[2], "--top") == 0) {
        for (const auto& prg_ptr : prgs) {
            npath = prg_ptr->prg.top_path();
            std::cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    } else if (strcmp(argv[2], "--bottom") == 0) {
        for (const auto& prg_ptr : prgs) {
            npath = prg_ptr->prg.bottom_path();
            std::cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    }

    // for each read in readfile,  infer node path along sequence
    std::string read_string;
    FastaqHandler readfile(argv[2]);
    while (not readfile.eof()) {
        try {
            readfile.get_next();
        } catch (std::out_of_range& err) {
            break;
        }

        for (const auto& prg_ptr : prgs) {
            npath = prg_ptr->prg.nodes_along_string(readfile.read);
            if (not npath.empty()) {
                std::cout << readfile.name << "\t" << prg_ptr->name << "\t";
                for (uint32_t j = 0; j != npath.size(); ++j) {
                    std::cout << "->" << npath[j]->id;
                }
                std::cout << std::endl;
            }
        }
    }

    return 0;
}
