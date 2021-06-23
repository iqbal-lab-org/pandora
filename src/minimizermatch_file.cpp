#include "minimizermatch_file.h"
#include "minihits.h"
#include "minihit.h"

MinimizerMatchFile::MinimizerMatchFile(const fs::path &filepath,
                                       const std::vector<std::shared_ptr<LocalPRG>>& prgs) :
    GenericFile(filepath), prgs(prgs) {
    file_handler << "read\tread_start\tread_end\tread_strand\tprg\tprg_mini_path" << std::endl;
}

void MinimizerMatchFile::write_hits(const Seq &seq, const Hits &hits) {
    std::vector<MinimizerHitPtr> sorted_hits;
    sorted_hits.insert(sorted_hits.end(), hits.begin(), hits.end());
    std::sort(sorted_hits.begin(), sorted_hits.end(), pCompReadPositionFirst());
    for (const MinimizerHitPtr &hit : hits) {
        file_handler << seq.name << "\t"
                     << hit->get_read_start_position() << "\t"
                     << hit->get_read_start_position()+hit->get_prg_path().length() << "\t"
                     << hit->read_strand << "\t"
                     << prgs[hit->get_prg_id()]->name << "\t"
                     << hit->get_prg_path() << std::endl;
    }
}
