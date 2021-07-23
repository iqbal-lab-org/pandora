#ifndef PANDORA_SAM_FILE_H
#define PANDORA_SAM_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"
#include "seq.h"
#include "localPRG.h"

class SAMFile : public GenericFile {
private:
    std::vector<bool> get_mapped_positions_bitset(const Seq &seq, const Hits &cluster) const;
    std::string get_cigar(const std::vector<bool> &mapped_positions_bitset) const;
    std::string get_segment_sequence(const Seq &seq,
                                     const std::vector<bool> &mapped_positions_bitset) const;
    const std::vector<std::shared_ptr<LocalPRG>>& prgs;
public:
    SAMFile(const fs::path &filepath,
            // just to convert prg IDs to prg names
            const std::vector<std::shared_ptr<LocalPRG>>& prgs);
    void write_sam_record_from_hit_cluster(
        const Seq &seq, const MinimizerHitClusters &clusters);
};

#endif // PANDORA_SAM_FILE_H
