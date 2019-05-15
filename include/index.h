#ifndef __INDEX_H_INCLUDED__   // if index.h hasn't been included yet...
#define __INDEX_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <memory>
#include "minirecord.h"
#include "prg/path.h"
#include "utils.h"


class Index {
public:
    std::unordered_map<uint64_t, std::vector<MiniRecord> *> minhash; //map of minimizers to MiniRecords - for each minimizer, records some information of it

    //declares all default constructors, destructors and assignment operators explicitly
    Index() = default; //default constructor
    Index(const Index &other) = default; //copy default constructor
    Index(Index &&other) = default; //move default constructor
    Index& operator=(const Index& other) = default; //copy assignment operator
    Index& operator=(Index&& other) = default; //move assignment operator
    virtual ~Index() = default; //destructor

    void add_record(const uint64_t, const uint32_t, const prg::Path, const uint32_t, const bool);

    void save(const std::string &prgfile, uint32_t w, uint32_t k);

    void save(const std::string &indexfile);

    void load(const std::string &prgfile, uint32_t w, uint32_t k);

    void load(const std::string &indexfile);

    void clear();

    bool operator==(const Index &other) const;

    bool operator!=(const Index &other) const;
};

void index_prgs(std::vector<std::shared_ptr<LocalPRG>> &,
                std::shared_ptr<Index> &,
                const uint32_t,
                const uint32_t,
                const std::string &,
                uint32_t threads=1);
#endif
