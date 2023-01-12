#include "localPRG_reader.h"

std::shared_ptr<LocalPRG> LocalPRGReaderGenerator::operator()() {
    while (!fh.eof()) {
        try {
            fh.get_next();
        } catch (std::out_of_range& err) {
            return nullptr;
        }
        if (fh.name.empty() or fh.read.empty())
            continue;
        auto local_prg = std::make_shared<LocalPRG>(LocalPRG(id, fh.name, fh.read));
        if (local_prg != nullptr) {
            id++;
            return local_prg;
        } else {
            fatal_error("Failed to make LocalPRG for ", fh.name);
        }
    }
    return nullptr;
}



std::vector<std::shared_ptr<LocalPRG>> read_prg_file_as_vector(const fs::path &filepath)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading PRGs from file " << filepath;
    std::vector<std::shared_ptr<LocalPRG>> prgs;

    FastaqHandler fh(filepath.string());
    uint32_t id = 0;
    while (!fh.eof()) {
        try {
            fh.get_next();
        } catch (std::out_of_range& err) {
            break;
        }
        if (fh.name.empty() or fh.read.empty())
            continue;
        auto s = std::make_shared<LocalPRG>(LocalPRG(id, fh.name,
            fh.read)); // build a node in the graph, which will represent a LocalPRG
                       // (the graph is a list of nodes, each representing a LocalPRG)
        if (s != nullptr) {
            prgs.push_back(s);
            id++;
        } else {
            fatal_error("Failed to make LocalPRG for ", fh.name);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Number of LocalPRGs read: " << prgs.size();

    return prgs;
}