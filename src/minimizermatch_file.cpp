#include "minimizermatch_file.h"
#include "minihits.h"
#include "minihit.h"
#include "denovo_discovery/denovo_utils.h"

MinimizerMatchFile::MinimizerMatchFile(const fs::path &filepath,
                                       const std::vector<std::string> &prg_names,
                                       bool is_fake_file) :
    GenericFile(filepath, is_fake_file), prg_names(prg_names) {
    (*this) << "kmer\tread\tread_start\tread_end\tread_strand\tprg\tprg_path\tprg_strand\n";
}

void MinimizerMatchFile::write_hits(const Seq &seq, const MinimizerHits &hits) {
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
        const std::string prg_name = prg_names[hit->get_prg_id()];
        (*this)  << kmer << "\t"
                 << seq.name << "\t"
                 << read_start_position << "\t"
                 << read_end_position << "\t"
                 << "-+"[hit->read_strand] << "\t"
                 << prg_name << "\t"
                 << hit->get_prg_path() << "\t"
                 << "-+"[hit->get_prg_kmer_strand()] << "\n";
    }
}
