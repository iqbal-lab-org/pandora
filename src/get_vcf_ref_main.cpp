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
#include "fastaq.h"


int pandora_get_vcf_ref(int argc, char *argv[]) // the "pandora walk" comand
{
    if (argc != 3 and argc != 2) {
        fprintf(stderr, "Usage: pandora get_vcf_ref <in_prg.fa> <seq.fa>\n");
        return 1;
    }

    // load prgs
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, argv[1]);

    // create fasta
    Fastaq fa(true, false);

    if (argc == 2) {
        for (const auto &prg_ptr: prgs) {
            std::vector<LocalNodePtr> npath;
            npath = prg_ptr->prg.top_path();
            fa.add_entry(prg_ptr->name, prg_ptr->string_along_path(npath));
        }
    } else {
        std::vector<LocalNodePtr> npath;
        std::string read_string;
        FastaqHandler readfile(argv[2]);
        bool found;

        for (const auto &prg_ptr: prgs) {
            found = false;
            readfile.get_id(0);
            while (not readfile.eof()) {
                npath = prg_ptr->get_valid_vcf_reference(readfile.read);
                if (not npath.empty()) {
                    //BOOST_LOG_TRIVIAL(debug) << ">" << prg_ptr->name << std::endl << readfile.read;
                    fa.add_entry(prg_ptr->name, prg_ptr->string_along_path(npath));
                    found = true;
                    break;
                }
                readfile.get_next();
            }

            if (!found) {
                assert(npath.empty());
                npath = prg_ptr->prg.top_path();
                fa.add_entry(prg_ptr->name, prg_ptr->string_along_path(npath));
            }
        }
    }

    std::string prg_file(argv[1]);
    fa.save(prg_file + ".vcf_ref.fa.gz");

    return 0;
}

