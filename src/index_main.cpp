#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>

#include "utils.h"
#include "localPRG.h"
#include <omp.h>

static void show_index_usage() {
    std::cerr << "Usage: pandora index [options] <prgs.fa>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              //<< "\t-u, --update\t\tLook for an index and add only PRGs with new names\n"
              << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, default 14\n"
              << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
              << "\t-t T\t\t\t\tNumber of threads, default 1\n"
              << "\t--offset\t\t\t\tOffset for PRG ids, default 0\n"
              << "\t--outfile\t\t\t\tFilename for index\n"
              << "\t--log_level\t\t\tdebug,[info],warning,error\n"
              << std::endl;
}

int pandora_index(int argc, char *argv[]) // the "pandora index" command
{
    // if not enough arguments, print usage
    if (argc < 2) {
        show_index_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    std::string prgfile, index_outfile, log_level = "info";
    bool update = false;
    uint32_t w = 14, k = 15, id=0, threads=1; // default parameters
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_index_usage();
            return 0;
        } else if ((arg == "-u") || (arg == "--update")) {
            update = true;
        } else if (arg == "-w") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                w = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-w option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-k") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                k = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-k option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-t") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                threads = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-t option requires one argument." << std::endl;
                return 1;
            }
        }
        else if (arg == "--offset") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                id = strtoul(argv[++i], nullptr, 10); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--offset option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "--outfile") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                index_outfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--outfile option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--log_level")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                log_level = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--log_level option requires one argument." << std::endl;
                return 1;
            }
        } else if (prgfile.empty()) {
            prgfile = argv[i];
            BOOST_LOG_TRIVIAL(debug) << "prgfile: " << prgfile;
        } else {
            std::cerr << argv[i] << " could not be attributed to any parameter" << std::endl;
        }
    }

    auto g_log_level{boost::log::trivial::info};
    if (log_level == "debug")
        g_log_level = boost::log::trivial::debug;
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= g_log_level);

    // load PRGs from file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, prgfile, id);

    // get output directory for the gfa
    boost::filesystem::path p(prgfile);
    boost::filesystem::path dir = p.parent_path();
    std::string outdir = dir.string();
    if (outdir.empty())
        outdir = ".";
    outdir += "/kmer_prgs";

    // index PRGs
    auto index = std::make_shared<Index>();
    index_prgs(prgs, index, w, k, outdir, threads);

    // save index
    BOOST_LOG_TRIVIAL(info) << "Saving index...";
    if (not index_outfile.empty())
        index->save(index_outfile);
    else if (id > 0)
        index->save(prgfile+"."+std::to_string(id), w, k);
    else
        index->save(prgfile, w, k);

    BOOST_LOG_TRIVIAL(info) << "All done!";
    return 0;
}

