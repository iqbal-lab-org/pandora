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

namespace fs = boost::filesystem;

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
public:
    std::unordered_map<uint64_t, std::vector<MiniRecord>*>
        minhash; // map of minimizers to MiniRecords - for each minimizer, records some
                 // information of it

    // declares all default constructors, destructors and assignment operators
    // explicitly
    Index() = default; // default constructor
    Index(const Index& other) = default; // copy default constructor
    Index(Index&& other) = default; // move default constructor
    Index& operator=(const Index& other) = default; // copy assignment operator
    Index& operator=(Index&& other) = default; // move assignment operator
    virtual ~Index() = default; // destructor

    void add_record(const AddRecordToIndexParams &params) {
        add_record(params.kmer, params.prg_id, params.path, params.knode_id,
            params.strand);
    }
    void add_record(
        const uint64_t, const uint32_t, const prg::Path&, const uint32_t, const bool);

    void save(const fs::path& prgfile, uint32_t w, uint32_t k);

    void save(const fs::path& indexfile);

    void load(fs::path prgfile, uint32_t w, uint32_t k);

    void load(const fs::path& indexfile);

    void clear();

    bool operator==(const Index& other) const;

    bool operator!=(const Index& other) const;

    void index_prgs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
        const uint32_t w, const uint32_t k, const fs::path& outdir,
        const uint32_t indexing_upper_bound=INDEXING_UPPER_BOUND_DEFAULT,
        uint32_t threads = 1);
};

#endif
