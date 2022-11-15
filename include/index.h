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
    const uint32_t w;
    const uint32_t k;
    ZipFile index_archive;
    std::vector<uint32_t> prg_min_path_lengths; // required to map to PRGs without actually loading them
    std::vector<uint32_t> prg_names;  // required to produce mapping files

    void save_minhash();

public:
    std::unordered_map<uint64_t, std::vector<MiniRecord>*> minhash; // TODO: move to private

    Index(const uint32_t w, const uint32_t k, const fs::path &outfilepath) :
        w(w), k(k), index_archive(outfilepath), minhash() {}
    virtual ~Index() {
        clear();
    };

    Index(Index&& other) = default; // move default constructor
    Index& operator=(Index&& other) = delete; // move assignment operator

    // not allowed to copy/assign - better to prohibit these due to RAM issue of
    // duplicating indexes
    Index(const Index& other) = delete;
    Index& operator=(const Index& other) = delete;

    void clear();

    void add_record(const AddRecordToIndexParams &params) {
        add_record(params.kmer, params.prg_id, params.path, params.knode_id,
            params.strand);
    }
    void add_record(
        const uint64_t, const uint32_t, const prg::Path&, const uint32_t, const bool);

    static Index load(const fs::path& indexfile);

    bool operator==(const Index& other) const;

    bool operator!=(const Index& other) const;

    /**
     * Index the given PRGs. Non-trivial parameters explanation given below.
     * @param optimise_RAM : a crucial parameter that allows us to optimise RAM when
     * indexing the PRGs, but has the critical side-effect of deleting the PRGs as we
     * scan through them. The index that is saved to the disk is
     * @param indexing_upper_bound
     */
    void index_prgs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        const uint32_t threads = 1);

    void index_prgs(LocalPRGReaderGeneratorIterator &prg_it,
        uintmax_t estimated_index_size=ESTIMATED_INDEX_SIZE_DEFAULT,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        const uint32_t threads=1);
};

#endif
