#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>

#include "utils.h"
#include "localPRG.h"




static void show_index_usage() {
    std::cerr << "Usage: pandora merge_index <index1> <index2> ...\n"
              << "Options:\n"
              << "\t--outfile\t\t\t\tFilename for merged index\n"
              << std::endl;
}

int pandora_merge_index(int argc, char *argv[]) // the "pandora merge_index" comand
{
    // if not enough arguments, print usage
    if (argc < 2) {
        show_index_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    std::string outfile = "merged_index.idx";
    std::vector<std::string> indexes;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_index_usage();
            return 0;
        } else if (arg == "--outfile") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                outfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else {// Uh-oh, there was no argument to the destination option.
                std::cerr << "--outfile option requires one argument." << std::endl;
                return 1;
            }
        } else {
            indexes.push_back(argv[i]);
        }
    }


    boost::filesystem::path p(outfile);
    boost::filesystem::path dir = p.parent_path();

    // merge indexes
    auto index = std::make_shared<Index>();
    for (const auto& new_index : indexes){
        index->load(new_index);
    }

    // save index
    index->save(outfile);

    return 0;
}

