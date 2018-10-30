#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>

#include "utils.h"
#include "localPRG.h"


using namespace std;

void index_prgs(std::vector<std::shared_ptr<LocalPRG>> &prgs, Index *idx, const uint32_t w, const uint32_t k,
                const string &outdir) {
    BOOST_LOG_TRIVIAL(debug) << "Index PRGs";

    // first reserve an estimated index size
    uint32_t r = 0;
    for (uint32_t i = 0; i != prgs.size(); ++i) {
        r += prgs[i]->seq.length();
    }
    idx->minhash.reserve(r);

    // now fill index
    auto dir_num = 0;
    for (uint32_t i = 0; i != prgs.size(); ++i) {
        if (i % 4000 == 0) {
            make_dir(outdir + "/" + int_to_string(dir_num + 1));
            dir_num++;
        }
        prgs[i]->minimizer_sketch(idx, w, k);
        prgs[i]->kmer_prg.save(
                outdir + "/" + int_to_string(dir_num) + "/" + prgs[i]->name + ".k" + to_string(k) + ".w" +
                to_string(w) + ".gfa");
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished adding " << prgs.size() << " LocalPRGs";
    BOOST_LOG_TRIVIAL(debug) << "Number of keys in Index: " << idx->minhash.size();
}

static void show_index_usage() {
    std::cerr << "Usage: pandora index [options] <prgs.fa>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              //<< "\t-u, --update\t\tLook for an index and add only PRGs with new names\n"
              << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, default 14\n"
              << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
              << std::endl;
}

int pandora_index(int argc, char *argv[]) // the "pandora index" comand
{
    // if not enough arguments, print usage
    if (argc < 2) {
        show_index_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile;
    bool update = false;
    uint32_t w = 14, k = 15; // default parameters
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_index_usage();
            return 0;
        } else if ((arg == "-u") || (arg == "--update")) {
            update = true;
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
        } else if (prgfile.empty()) {
            prgfile = argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
            BOOST_LOG_TRIVIAL(debug) << "prgfile: " << prgfile;
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }


    // load PRGs from file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, prgfile);

    // get output directory for the gfa
    boost::filesystem::path p(prgfile);
    boost::filesystem::path dir = p.parent_path();
    string outdir = dir.string();
    if (outdir.empty())
        outdir = ".";
    outdir += "/kmer_prgs";

    // index PRGs
    Index *idx;
    idx = new Index();
    index_prgs(prgs, idx, w, k, outdir);

    // save index
    idx->save(prgfile, w, k);

    return 0;
}

