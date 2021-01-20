#ifndef __UTILS_H_INCLUDED__ // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

#include <vector>
#include <set>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <limits>
#include <boost/filesystem/path.hpp>
#include "minihits.h"
#include "pangenome/ns.cpp"
#include <boost/log/trivial.hpp>
#include <sstream>
#include "fatal_error.h"

namespace fs = boost::filesystem;

namespace pandora {
enum class Strand {
    Forward = '+',
    Reverse = '-',
};
}
class Index;

class PanNode;

class LocalPRG;

class Seq;

typedef std::unordered_map<std::string, std::string> VCFRefs;

template <typename T> struct pointer_values_equal {
    const T* to_find;

    bool operator()(const T* other) const { return *to_find == *other; }
};

template <typename T> struct spointer_values_equal {
    const std::shared_ptr<T> to_find;

    bool operator()(const std::shared_ptr<T> other) const { return *to_find == *other; }
};

// pointer less comparator, from
// https://stackoverflow.com/questions/41375232/is-there-an-stl-comparator-for-stdset-or-stdmap-with-shared-ptr-keys-that
struct ptr_less {
    template <typename T> bool operator()(T lhs, T rhs) const
    {
        return std::less<decltype(*lhs)>()(*lhs, *rhs);
    }
};

// utility functions
std::string now();

std::string int_to_string(const int number);

std::vector<std::string> split(const std::string&, const std::string&);

char complement(char);

std::string rev_complement(std::string);

float lognchoosek2(uint32_t, uint32_t, uint32_t);

// probably should be moved to map_main.cpp
void read_prg_file(std::vector<std::shared_ptr<LocalPRG>>& prgs,
    const fs::path& filepath, uint32_t id = 0);

void load_PRG_kmergraphs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
    const uint32_t& w, const uint32_t& k, const fs::path& prgfile);

void load_vcf_refs_file(const fs::path& filepath, VCFRefs& vcf_refs);

void add_read_hits(const Seq&, const std::shared_ptr<MinimizerHits>&, const Index&);

void define_clusters(std::set<std::set<MinimizerHitPtr, pComp>, clusterComp>&,
    const std::vector<std::shared_ptr<LocalPRG>>&, std::shared_ptr<MinimizerHits>,
    const int, const float&, const uint32_t, const uint32_t);

void filter_clusters(std::set<std::set<MinimizerHitPtr, pComp>, clusterComp>&);

void filter_clusters2(
    std::set<std::set<MinimizerHitPtr, pComp>, clusterComp>&, const uint32_t&);

void infer_localPRG_order_for_reads(const std::vector<std::shared_ptr<LocalPRG>>& prgs,
    std::shared_ptr<MinimizerHits> minimizer_hits, std::shared_ptr<pangenome::Graph>,
    const int, const uint32_t&, const float&, const uint32_t min_cluster_size = 10,
    const uint32_t expected_number_kmers_in_short_read_sketch
    = std::numeric_limits<uint32_t>::max());

uint32_t pangraph_from_read_file(const std::string&, std::shared_ptr<pangenome::Graph>,
    std::shared_ptr<Index>, const std::vector<std::shared_ptr<LocalPRG>>&,
    const uint32_t, const uint32_t, const int, const float&,
    const uint32_t min_cluster_size = 10, const uint32_t genome_size = 5000000,
    const bool illumina = false, const bool clean = false,
    const uint32_t max_covg = 300, uint32_t threads = 1);

//, const uint32_t, const float&, bool);
void infer_most_likely_prg_path_for_pannode(
    const std::vector<std::shared_ptr<LocalPRG>>&, PanNode*, uint32_t, float);

// TODO : refactor all file open and closing to use these functions
void open_file_for_reading(const std::string& file_path, std::ifstream& stream);
void open_file_for_writing(const std::string& file_path, std::ofstream& stream);

std::vector<std::string> get_vector_of_strings_from_file(const std::string& file_path);

// string to genome size
// effectively a copy of
// https://github.com/lh3/minimap2/blob/6a4b9f9082b66597185a97c847b548250363d65a/main.c#L84
uint32_t strtogs(const char*);

// used to transform the CLI string
// https://cliutils.github.io/CLI11/book/chapters/validators.html
std::string transform_cli_gsize(std::string);

// used to transform paths to absolute paths - designed to be used with CLI11 transform
std::string make_absolute(std::string);

#endif
