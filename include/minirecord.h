#ifndef __MINIRECORD_H_INCLUDED__   // if minirecord.h hasn't been included yet...
#define __MINIRECORD_H_INCLUDED__

#include <iostream>
#include <cstdint>
#include "prg/path.h"

//Minimizer Record
struct MiniRecord {
    uint32_t prg_id; //prg id of the minimizer
    prg::Path path; //kmer path of the minimizer
    uint32_t knode_id; //kmer graph node id
    bool strand;

    MiniRecord();

    MiniRecord(const uint32_t, const prg::Path, const uint32_t, const bool);

    ~MiniRecord();

    bool operator==(const MiniRecord &y) const;

    friend std::ostream &operator<<(std::ostream &out, const MiniRecord &m);

    friend std::istream &operator>>(std::istream &in, MiniRecord &m);
};

#endif
