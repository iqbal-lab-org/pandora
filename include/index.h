#ifndef __INDEX_H_INCLUDED__   // if index.h hasn't been included yet...
#define __INDEX_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include "minirecord.h"
#include "prg/path.h"


class Index {
public:
    std::unordered_map<uint64_t, std::vector<MiniRecord> *> minhash;

    Index();

    ~Index();

    void add_record(const uint64_t, const uint32_t, const prg::Path, const uint32_t, const bool);

    void save(const std::string &prgfile, uint32_t w, uint32_t k);

    void load(const std::string &prgfile, uint32_t w, uint32_t k);

    void clear();
};

#endif
