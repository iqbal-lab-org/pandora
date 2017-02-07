#ifndef __INDEX_H_INCLUDED__   // if index.h hasn't been included yet...
#define __INDEX_H_INCLUDED__

#include <stdint.h>
#include <ostream>
#include <vector>
#include <map>
#include "minirecord.h"
#include "path.h"

class Index {
  public:
    std::map<uint64_t,std::vector<MiniRecord>> minhash;
    Index();
    ~Index();
    void add_record(uint64_t, uint32_t, Path, bool);
    void save_index(const std::string& filename);
    void clear();
  friend std::ostream& operator<< (std::ostream& out, const Index& idx);
};

#endif
