#include "paf_file.h"
#include "minihits.h"
#include "minihit.h"
PafFile::PafFile(const fs::path &filepath) : GenericFile(filepath) {
    file_handler << "read_idx\tmlen\tmstart\tmend\tmstrands\tprg_idx\tnmatch\talen\tmapq" << std::endl;
}

void PafFile::write_clusters(const MinimizerHitClusters &clusters) {
    uint32_t number_of_hits_in_read = 0;
    for (const Hits &cluster : clusters) {
        number_of_hits_in_read += cluster.size();
    }

    std::vector<MinimizerHitPtr> sorted_minimizer_hits;
    for (const Hits&cluster : clusters) {
        sorted_minimizer_hits.insert(sorted_minimizer_hits.end(),
            cluster.begin(), cluster.end());
    }
    std::sort(sorted_minimizer_hits.begin(), sorted_minimizer_hits.end(), pComp());

    for (const Hits&cluster : clusters) {
        this->write_cluster(cluster, number_of_hits_in_read, sorted_minimizer_hits);
    }
}

void PafFile::write_cluster(const Hits &cluster, uint32_t number_of_hits_in_read,
                   const std::vector<MinimizerHitPtr> &sorted_minimizer_hits) {
    const MinimizerHitPtr first_hit = *(cluster.begin());
    std::vector<size_t> positions_of_hits;
    for (const MinimizerHitPtr &hit : cluster) {
        size_t position_of_this_hit =
            find(sorted_minimizer_hits.begin(), sorted_minimizer_hits.end(), hit) - sorted_minimizer_hits.begin();
        positions_of_hits.push_back(position_of_this_hit);
    }
    size_t lowest_position = *std::min_element(positions_of_hits.begin(), positions_of_hits.end());
    size_t highest_position = *std::max_element(positions_of_hits.begin(), positions_of_hits.end());

    std::vector<char> strands;
    for (const MinimizerHitPtr &hit : cluster) {
        if (hit->read_strand)
            strands.push_back('+');
        else
            strands.push_back('-');
    }
    std::string strands_info =
        to_string(std::count(strands.begin(), strands.end(), '+')) + "+/" +
        to_string(std::count(strands.begin(), strands.end(), '-')) + "-";

    size_t nmatch = cluster.size();
    size_t alen = highest_position - lowest_position + 1;
    uint32_t mapq = 256 - (nmatch / alen * 256);

    file_handler << first_hit->get_read_id() << "\t"
                 << number_of_hits_in_read << "\t"
                 << lowest_position << "\t"
                 << highest_position << "\t"
                 << strands_info << "\t"
                 << first_hit->minimizer_from_PRG.prg_id << "\t"
                 << nmatch << "\t"
                 << alen << "\t"
                 << mapq << std::endl;
}
