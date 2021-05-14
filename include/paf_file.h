#ifndef PANDORA_PAF_FILE_H
#define PANDORA_PAF_FILE_H

#include "forward_declarations.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "utils.h"
#include <fstream>
#include "minihit.h"
#include <algorithm>

class PafFile {
private:
    std::ofstream paf_file_handler;
    void write_cluster(const MinimizerHitCluster &cluster, uint32_t number_of_hits_in_read,
        const std::vector<MinimizerHitPtr> &sorted_minimizer_hits);

public:
    PafFile(const fs::path &filepath);
    virtual ~PafFile();
    void write_clusters(const MinimizerHitClusters &clusters);
    void write_cluster(const MinimizerHitCluster &clusters);
};

#endif // PANDORA_PAF_FILE_H
