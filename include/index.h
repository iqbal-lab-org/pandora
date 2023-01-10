#ifndef __INDEX_H_INCLUDED__ // if index.h hasn't been included yet...
#define __INDEX_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>
#include <boost/filesystem.hpp>
#include "minirecord.h"
#include "prg/path.h"
#include "utils.h"
#include "globals.h"
#include "zip_file.h"
#include "localPRG_reader.h"
#include "forward_declarations.h"

struct AddRecordToIndexParams {
    const uint64_t kmer;
    const uint32_t prg_id;
    const prg::Path path;
    const uint32_t knode_id;
    const bool strand;
    AddRecordToIndexParams(const uint64_t kmer, const uint32_t prg_id,
        const prg::Path& path, const uint32_t knode_id, const bool strand) :
        kmer(kmer), prg_id(prg_id), path(path), knode_id(knode_id), strand(strand){ }
};

class Index {
private:
    ///////////////////////////////////////////////////////////////////////////////////
    // attributes
    const uint32_t w;
    const uint32_t k;

    // required to produce mapping files
    std::vector<std::string> prg_names;

    // required to map to PRGs without actually loading them
    std::vector<uint32_t> prg_lengths;
    std::vector<uint32_t> prg_min_path_lengths;

    // populated when loading an index - just valid for index consumption commands
    // e.g. map, compare, discover, etc
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    ///////////////////////////////////////////////////////////////////////////////////

    // By making constructors private we are explicitly disallowing the index to be
    // built from other classes/modules
    // Note: prg_names and prg_lengths are r-value refs because we take ownership of it when building an index
    Index(const uint32_t w, const uint32_t k, std::vector<std::string> &&prg_names,
        std::vector<uint32_t> &&prg_lengths) :
        w(w), k(k), prg_names(std::move(prg_names)), prg_lengths(std::move(prg_lengths)) {
        prg_min_path_lengths.reserve(get_number_of_prgs());
        prgs.reserve(get_number_of_prgs());
    }
    Index(const uint32_t w, const uint32_t k, std::vector<std::string> &&prg_names,
        std::vector<uint32_t> &&prg_lengths, std::vector<uint32_t> &&prg_min_path_lengths) :
        w(w), k(k), prg_names(std::move(prg_names)), prg_lengths(std::move(prg_lengths)),
        prg_min_path_lengths(std::move(prg_min_path_lengths)) {
        prgs.reserve(get_number_of_prgs());
    }

    void load_prgs(ZipFileReader &zip_file);

    void index_prgs(ZipFileWriter &index_archive,
        LocalPRGReaderGeneratorIterator &prg_it,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        const uint32_t threads=1);

    std::unordered_map<std::string, uint32_t> get_prg_names_to_prg_index() const;

    void save_minhash(ZipFileWriter &index_archive) const;
    void load_minhash(ZipFileReader &zip_file);

    template <class Iterator>
    static void save_values(ZipFileWriter &index_archive, const std::string &zip_path,
        Iterator begin, Iterator end) {
        index_archive.prepare_new_entry(zip_path);
        for (; begin != end; ++begin) {
            std::stringstream ss;
            ss << *begin << std::endl;
            index_archive.write_data(ss.str());
        }
    }
    inline void save_prg_names(ZipFileWriter &index_archive) const {
        save_values(index_archive, "_prg_names", prg_names.begin(), prg_names.end());
    }
    inline void save_prg_lengths(ZipFileWriter &index_archive) const {
        save_values(index_archive, "_prg_lengths",
            prg_lengths.begin(), prg_lengths.end());
    }
    inline void save_prg_min_path_lengths(ZipFileWriter &index_archive) const {
        save_values(index_archive, "_prg_min_path_lengths",
            prg_min_path_lengths.begin(), prg_min_path_lengths.end());
    }
    inline void save_metadata(ZipFileWriter &index_archive) const {
        std::vector<uint32_t> metadata{w, k};
        save_values(index_archive, "_metadata", metadata.begin(), metadata.end());
    }

public:
    ///////////////////////////////////////////////////////////////////////////////////
    // attributes
    std::unordered_map<uint64_t, std::vector<MiniRecord>*> minhash; // TODO: move to private
    ///////////////////////////////////////////////////////////////////////////////////

    virtual ~Index() {
        clear();
    };

    // moving indexes is allowed
    Index(Index&& other) = default;
    Index& operator=(Index&& other) = default;

    // copying indexes is not allowed
    Index(const Index& other) = delete;
    Index& operator=(const Index& other) = delete;

    ////////////////////////////////////////////////////////////////////////////////////
    // public interface to build and load indexes
    static void build_index_on_disk(const uint32_t w, const uint32_t k,
        const fs::path &prg_filepath, const fs::path &out_filepath,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        const uint32_t threads=1);

    /**
     * Load the index from a given path
     * @param indexfile : the path to the index file
     * @param lazy_load : if true, does not preemptively load all PRGs and kmer PRGs
     * @return the index
     */
    static Index load(const fs::path& indexfile, bool lazy_load = false);
    ////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////
    // misc getters
    inline uint32_t get_window_size() const {
        return w;
    }

    inline uint32_t get_kmer_size() const {
        return k;
    }

    inline size_t get_number_of_prgs() const {
        return prg_names.size();
    }

    inline const std::vector<std::string>& get_prg_names() const {
        return prg_names;
    }
    inline const std::string & get_prg_name_given_id(size_t id) const {
        return prg_names[id];
    }

    inline const std::vector<uint32_t>& get_prg_lengths() const {
        return prg_lengths;
    }
    inline uint32_t get_prg_length_given_id(size_t id) const {
        return prg_lengths[id];
    }

    inline const std::vector<std::shared_ptr<LocalPRG>> & get_prgs() const {
        return prgs;
    }
    inline const std::shared_ptr<LocalPRG>& get_prg_given_id(size_t id) const {
        return prgs[id];
    }
    ////////////////////////////////////////////////////////////////////////////////////

    void clear();

    void add_record(const AddRecordToIndexParams &params) {
        add_record(params.kmer, params.prg_id, params.path, params.knode_id,
            params.strand);
    }
    void add_record(
        const uint64_t, const uint32_t, const prg::Path&, const uint32_t, const bool);

    bool operator==(const Index& other) const;

    bool operator!=(const Index& other) const;
};

#endif
