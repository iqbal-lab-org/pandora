#include "fastaq.h"
#include "fastaq_handler.h"
#include "localPRG.h"
#include "localgraph.h"
#include "localnode.h"
#include "utils.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

int pandora_random_path(int argc, char* argv[]) // the "pandora walk" comand
{
    if (argc != 3 and argc != 2) {
        fprintf(stderr, "Usage: pandora random_path <in_prg.fa> [<num_paths>]\n");
        return 1;
    }

    // load prgs
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, argv[1]);

    // create fasta
    Fastaq fa(true, false);

    uint32_t num_paths = 1;
    if (argc == 3) {
        num_paths = strtoul(argv[2], nullptr, 10);
    }

    for (const auto& prg_ptr : prgs) {
        std::unordered_set<std::string> random_paths;
        auto skip = 0;
        while (random_paths.size() < num_paths and skip < 10) {
            auto spath = prg_ptr->random_path();
            if (random_paths.find(spath) != random_paths.end()) {
                skip += 1;
            } else {
                random_paths.insert(spath);
            }
        }
        uint32_t i = 0;
        for (const auto& path : random_paths) {
            fa.add_entry(prg_ptr->name + "_" + std::to_string(i), path);
            i++;
        }
    }

    fa.save("random_paths.fa.gz");

    return 0;
}
