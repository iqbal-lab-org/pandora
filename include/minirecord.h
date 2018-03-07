#ifndef __MINIRECORD_H_INCLUDED__   // if minirecord.h hasn't been included yet...
#define __MINIRECORD_H_INCLUDED__

#include <iostream>
#include <algorithm>
#include "path.h"

struct MiniRecord {
    uint32_t prg_id;
    Path path;
    uint32_t knode_id;
    bool strand;

    MiniRecord();

    MiniRecord(const uint32_t, const Path, const uint32_t, const bool);

    ~MiniRecord();

    bool operator==(const MiniRecord &y) const;

    friend std::ostream &operator<<(std::ostream &out, const MiniRecord &m);

    friend std::istream &operator>>(std::istream &in, MiniRecord &m);
};

#endif
