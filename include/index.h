#ifndef __INDEX_H_INCLUDED__   // if index.h hasn't been included yet...
#define __INDEX_H_INCLUDED__

#include <stdint.h>
#include <ostream>
#include <vector>
#include <map>
#include "minirecord.h"
#include "path.h"

using namespace std;

class Index {
  public:
    map<uint64_t,vector<MiniRecord>> minhash;
    Index();
    ~Index();
    void add_record(uint64_t, uint32_t, Path);
    void save_index(string filename);
    void clear();
  friend ostream& operator<< (ostream& out, const Index& idx);
};

#endif
