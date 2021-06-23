#include "minimizermatch_file.h"
#include "minihits.h"
#include "minihit.h"

MinimizerMatchFile::MinimizerMatchFile(const fs::path &filepath) : GenericFile(filepath) {
    file_handler << "read_idx\tread_start\tread_end\tread_strand\tprg_idx\tprg_mini_path" << std::endl;
}

void MinimizerMatchFile::write_hits(const Hits &hits) {
    std::vector<MinimizerHitPtr> sorted_hits;
    sorted_hits.insert(sorted_hits.end(), hits.begin(), hits.end());
    std::sort(sorted_hits.begin(), sorted_hits.end(), pCompReadPositionFirst());
    for (const MinimizerHitPtr &hit : hits) {
        file_handler << hit->get_read_id() << "\t"
                     << hit->get_read_start_position() << "\t"
                     << hit->get_read_start_position()+hit->get_prg_path().length() << "\t"
                     << hit->read_strand << "\t"
                     << hit->get_prg_id() << "\t"
                     << hit->get_prg_path() << std::endl;
    }
}
