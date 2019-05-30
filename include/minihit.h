#ifndef __MINIHIT_H_INCLUDED__   // if minihit.h hasn't been included yet...
#define __MINIHIT_H_INCLUDED__

#include <ostream>
#include <cstdint>
#include "minimizer.h"
#include "minirecord.h"

//TODO: here we have one MinimizerHit for each (read_id, read_start_position, read_strand) and MiniRecord
//TODO: we could make one (read_id, read_start_position, read_strand) and a vector of MiniRecord
struct MinimizerHit {
private:
    uint32_t read_id; //TODO: this can be made a template and change depending on the maximum number of reads
    uint32_t read_start_position; //TODO: this can be made a template and change depending on the maximum read length
    bool read_strand;
    const MiniRecord &minimizerFromPRG;

public:
    inline uint32_t get_read_id() const { return read_id; }
    inline uint32_t get_read_start_position() const  { return read_start_position; }
    inline uint32_t get_prg_id () const { return minimizerFromPRG.prg_id; }
    inline const prg::Path & get_prg_path() const { return minimizerFromPRG.path; }
    inline uint32_t get_kmer_node_id() const { return minimizerFromPRG.knode_id; }
    inline bool is_forward() const { return read_strand == minimizerFromPRG.strand ;}

    MinimizerHit(const uint32_t i, const Minimizer &minimizerFromRead, const MiniRecord &minimizerFromPRG);

    bool operator<(const MinimizerHit &y) const;

    bool operator==(const MinimizerHit &y) const;

    friend std::ostream &operator<<(std::ostream &out, const MinimizerHit &m);
};


#endif
