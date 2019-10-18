#ifndef __MINIHIT_H_INCLUDED__ // if minihit.h hasn't been included yet...
#define __MINIHIT_H_INCLUDED__

#include "minimizer.h"
#include "minirecord.h"
#include <cstdint>
#include <ostream>

/**
 * Describes a hit between a read an a minimizer from the PRG
 * TODO: Possible improvement (memory): here we have one MinimizerHit for each (read_id,
 * read_start_position, read_strand) and MiniRecord
 * TODO: Possible improvement (memory): we could make one (read_id, read_start_position,
 * read_strand) and a vector of MiniRecords
 */
struct MinimizerHit {
private:
    uint32_t read_id; // TODO: Possible improvement (memory): this can be made a
                      // template and change depending on the maximum number of reads
    uint32_t
        read_start_position; // TODO: Possible improvement (memory): this can be made a
                             // template and change depending on the maximum read length
    bool read_strand;
    const MiniRecord& minimizer_from_PRG;

public:
    inline uint32_t get_read_id() const { return read_id; }
    inline uint32_t get_read_start_position() const { return read_start_position; }
    inline uint32_t get_prg_id() const { return minimizer_from_PRG.prg_id; }
    inline const prg::Path& get_prg_path() const { return minimizer_from_PRG.path; }
    inline uint32_t get_kmer_node_id() const { return minimizer_from_PRG.knode_id; }
    inline bool is_forward() const
    {
        return read_strand == minimizer_from_PRG.strand;
    } // TODO: the name of this method is very misleading, should be same_strands() or
      // sth like this

    MinimizerHit(const uint32_t i, const Minimizer& minimizer_from_read,
        const MiniRecord& minimizer_from_PRG);

    bool operator<(const MinimizerHit& y) const;

    bool operator==(const MinimizerHit& y) const;

    friend std::ostream& operator<<(std::ostream& out, const MinimizerHit& m);
};

#endif
