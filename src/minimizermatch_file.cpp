#include "minimizermatch_file.h"
#include "minihits.h"
#include "minihit.h"
#include "denovo_discovery/denovo_utils.h"

MinimizerMatchFile::MinimizerMatchFile(const fs::path &filepath,
                                       const std::vector<std::shared_ptr<LocalPRG>>& prgs,
                                       bool is_fake_file) :
    GenericFile(filepath, is_fake_file), prgs(prgs) {
    (*this) << "kmer\tread\tread_start\tread_end\tread_strand\tprg\tprg_strand\tprg_path\n";
}

void MinimizerMatchFile::write_hits(const Seq &seq, const Hits &hits) {
    std::vector<MinimizerHitPtr> sorted_hits;
    sorted_hits.insert(sorted_hits.end(), hits.begin(), hits.end());
    std::sort(sorted_hits.begin(), sorted_hits.end(), pCompReadPositionFirst());
    for (const MinimizerHitPtr &hit : hits) {
        const uint32_t read_start_position = hit->get_read_start_position();
        const uint32_t read_end_position = hit->get_read_start_position()+hit->get_prg_path().length();
        std::string kmer = seq.full_seq.substr(read_start_position, hit->get_prg_path().length());
        if (hit->read_strand == 0) {
            kmer = reverse_complement(kmer);
        }
        const std::string prg_name = prgs[hit->get_prg_id()]->name;
        (*this)  << kmer << "\t"
                 << seq.name << "\t"
                 << read_start_position << "\t"
                 << read_end_position << "\t"
                 << hit->read_strand << "\t"
                 << prg_name << "\t"
                 << hit->get_prg_kmer_strand() << "\t"
                 << hit->get_prg_path() << "\n";
    }
}
