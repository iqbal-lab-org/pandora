#ifndef __MINIHIT_H_INCLUDED__ // if minihit.h hasn't been included yet...
#define __MINIHIT_H_INCLUDED__

#include <ostream>
#include <cstdint>
#include "minimizer.h"
#include "minirecord.h"
#include "fatal_error.h"

/**
 * Describes a hit between a read an a minimizer from the PRG
 * TODO: Possible improvement (memory): here we have one MinimizerHit for each (read_id,
 * read_start_position, read_strand) and MiniRecord
 * TODO: Possible improvement (memory): we could make one (read_id, read_start_position,
 * read_strand) and a vector of MiniRecords
 */
struct MinimizerHit {
private:
    uint32_t read_id;
    uint32_t read_start_position;
public:
    const bool read_strand;
    const MiniRecord& minimizer_from_PRG;
    inline uint32_t get_read_id() const { return read_id; }
    inline uint32_t get_read_start_position() const { return read_start_position; }
    inline uint32_t get_prg_id() const { return minimizer_from_PRG.prg_id; }
    inline const prg::Path& get_prg_path() const { return minimizer_from_PRG.path; }
    inline uint32_t get_kmer_node_id() const { return minimizer_from_PRG.knode_id; }
    inline bool get_prg_kmer_strand() const { return minimizer_from_PRG.strand; }
    inline bool same_strands() const
    {
        return read_strand == get_prg_kmer_strand();
    }

    MinimizerHit(const uint32_t i, const Minimizer& minimizer_from_read,
        const MiniRecord& minimizer_from_PRG);

    bool operator<(const MinimizerHit& y) const;

    bool operator==(const MinimizerHit& y) const;

    friend std::ostream& operator<<(std::ostream& out, const MinimizerHit& m);
};

#endif
