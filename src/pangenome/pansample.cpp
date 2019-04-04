#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <unordered_set>
#include <algorithm>
#include "pangenome/pansample.h"
#include "pangenome/pannode.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

Sample::Sample(const std::string &s, const uint32_t &id) : name(s), sample_id(id) {}

void Sample::add_path(const uint32_t node_id, const std::vector<KmerNodePtr> &c) {
    if (paths.find(node_id) == paths.end()) {
        paths[node_id] = {c};
    } else {
        paths[node_id].push_back(c);
    }
}

bool Sample::operator==(const Sample &y) const {
    if (name != y.name) { return false; }
    return true;
}

bool Sample::operator!=(const Sample &y) const {
    return !(*this == y);
}

bool Sample::operator<(const Sample &y) const {
    return (name < y.name);
}

std::ostream &pangenome::operator<<(std::ostream &out, pangenome::Sample const &s) {
    out << s.name << ":\t";
    for (const auto &p : s.paths) {
        for (uint32_t i = 0; i != p.second.size(); ++i) {
            out << p.first << "\t";
            for (const auto &n : p.second[i]){
                out << *n;
            }
            out << std::endl;
        }
    }
    return out;
}
