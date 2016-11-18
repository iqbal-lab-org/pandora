#ifndef __MINIRECORD_H_INCLUDED__   // if minirecord.h hasn't been included yet...
#define __MINIRECORD_H_INCLUDED__

#include <ostream>
#include <algorithm>
#include "path.h"

using namespace std;

struct MiniRecord
{
    uint32_t prg_id;
    Path path;
    bool strand;
    MiniRecord(const uint32_t, const Path, const bool);
    ~MiniRecord();
    bool operator == (const MiniRecord& y) const;
    friend ostream& operator<< (ostream& out, const MiniRecord& m);
};

#endif
