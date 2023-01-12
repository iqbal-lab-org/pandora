#ifndef __UTILS_H_INCLUDED__ // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

#include "forward_declarations.h"
#include <vector>
#include <set>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <string>
#include <limits>
#include <utility>
#include <boost/filesystem/path.hpp>
#include "minihits.h"
#include <boost/log/trivial.hpp>
#include <sstream>
#include "fatal_error.h"
#include "paf_file.h"
#include "inthash.h"
#include "cluster_files.h"

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

void load_vcf_refs_file(const fs::path& filepath, VCFRefs& vcf_refs);

void add_read_hits(const Seq&, const std::shared_ptr<MinimizerHits>&, const Index&);

void define_clusters(
    const std::string &sample_name,
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits,
    const std::vector<uint32_t> &prg_min_path_lengths,
    const std::vector<std::string> &prg_names,
    std::shared_ptr<MinimizerHits> minimizer_hits, const int max_diff,
    const float& fraction_kmers_required_for_cluster, const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch,
    ClusterDefFile& cluster_def_file);

void filter_clusters(
    const std::string &sample_name,
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits,
    const std::vector<std::string> &prg_names,
    ClusterFilterFile& cluster_filter_file
);

void add_clusters_to_pangraph(
    const MinimizerHitClusters& minimizer_hit_clusters,
    std::shared_ptr<pangenome::Graph> &pangraph,
    Index &index);

MinimizerHitClusters get_minimizer_hit_clusters(
    const std::string &sample_name,
    const Seq &seq,
    const std::vector<uint32_t> &prg_min_path_lengths,
    const std::vector<std::string> &prg_names,
    std::shared_ptr<MinimizerHits> minimizer_hits,
    std::shared_ptr<pangenome::Graph> pangraph, const int max_diff,
    const uint32_t& genome_size, const float& fraction_kmers_required_for_cluster,
    ClusterDefFile &cluster_def_file,
    ClusterFilterFile &cluster_filter_file,
    const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch = std::numeric_limits<uint32_t>::max());

uint32_t pangraph_from_read_file(const SampleData& sample,
    std::shared_ptr<pangenome::Graph> pangraph, Index &index,
    const int max_diff, const float& e_rate,
    const fs::path& sample_outdir, const uint32_t min_cluster_size = 10,
    const uint32_t genome_size = 5000000, const bool illumina = false, const bool clean = false,
    const uint32_t max_covg = 300, uint32_t threads = 1,
    const bool keep_extra_debugging_files = false);

void infer_most_likely_prg_path_for_pannode(
    const std::vector<std::shared_ptr<LocalPRG>>&, PanNode*, uint32_t, float);

// TODO : refactor all file open and closing to use these functions
void open_file_for_reading(const std::string& file_path, std::ifstream& stream);
void open_file_for_writing(const std::string& file_path, std::ofstream& stream);
void open_file_for_appending(const std::string& file_path, std::ofstream& stream);

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

std::vector<SampleData> load_read_index(const fs::path& read_index_fpath);

// Builds a file in memory
// Returns filepath
std::pair<int, std::string> build_memfd(const std::string &data);

std::string exec(const char* cmd);

void build_file(const std::string &filepath, const std::string &data);

bool tool_exists(const std::string &command);

void concatenate_text_files(
    const fs::path& output_filename, const std::vector<fs::path>& input_filenames,
    const std::string &prepend="");

std::string reverse_complement(const std::string& forward);

inline void to_upper(std::string &str) {
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

int random_int();
std::pair<std::vector<std::string>, std::vector<size_t>> split_ambiguous(const std::string& input_string, uint8_t delim = 4);

uintmax_t get_number_of_bytes_in_file(const fs::path &file);

#endif
