#include "paf_file.h"
#include "minihits.h"
#include "minihit.h"
#include "localPRG.h"

PafFile::PafFile(const fs::path &filepath,
                 const std::vector<std::shared_ptr<LocalPRG>>& prgs,
                 bool is_fake_file)
    : GenericFile(filepath, is_fake_file), prgs(prgs) {
    (*this) << "qname\tqmlen\tqmstart\tqmend\tqmstrands\tprg\tnmatch\talen\tmapq\n";
}

void PafFile::write_clusters(const Seq &seq, const MinimizerHitClusters &clusters) {
    std::vector<MinimizerHitPtr> all_sorted_minimizer_hits;
    for (const MinimizerHits& cluster : clusters) {
        all_sorted_minimizer_hits.insert(all_sorted_minimizer_hits.end(),cluster.begin(), cluster.end());
    }
    std::sort(
        all_sorted_minimizer_hits.begin(), all_sorted_minimizer_hits.end(), pComp());

    for (const MinimizerHits&cluster : clusters) {
        this->write_cluster(seq, cluster, all_sorted_minimizer_hits);
    }
}

void PafFile::write_cluster(const Seq &seq, const MinimizerHits &cluster,
                   const std::vector<MinimizerHitPtr> &all_sorted_minimizer_hits) {
    const MinimizerHitPtr first_hit = *(cluster.begin());
    std::vector<size_t> positions_of_hits;
    for (const MinimizerHitPtr &hit : cluster) {
        size_t position_of_this_hit =
            find(all_sorted_minimizer_hits.begin(), all_sorted_minimizer_hits.end(), hit) - all_sorted_minimizer_hits.begin();
        positions_of_hits.push_back(position_of_this_hit);
    }
    size_t lowest_position = *std::min_element(positions_of_hits.begin(), positions_of_hits.end());
    size_t highest_position = *std::max_element(positions_of_hits.begin(), positions_of_hits.end());

    std::vector<char> strands;
    for (const MinimizerHitPtr &hit : cluster) {
        strands.push_back("-+"[hit->same_strands()]);
    }
    uint32_t plus_strand_count = std::count(strands.begin(), strands.end(), '+');
    uint32_t minus_strand_count = std::count(strands.begin(), strands.end(), '-');

    std::string strands_info;
    if (minus_strand_count == 0) {
        // we are sure strand is +
        strands_info = "+";
    } else if (plus_strand_count == 0) {
        // we are sure strand is -
        strands_info = "-";
    } else {
        // we have some hits saying strand is + and other saying strand is -
        strands_info = "?(" +
            to_string(plus_strand_count) + "+/"  +
            to_string(minus_strand_count) + "-)";
    }

    size_t nmatch = cluster.size();
    size_t alen = highest_position - lowest_position + 1;
    uint32_t mapq = 255;

    (*this)  << seq.name << "\t"
             << all_sorted_minimizer_hits.size() << "\t"
             << lowest_position << "\t"
             << highest_position << "\t"
             << strands_info << "\t"
             << prgs[first_hit->get_prg_id()]->name << "\t"
             << nmatch << "\t"
             << alen << "\t"
             << mapq << "\n";
}
