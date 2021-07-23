#ifndef PANDORA_PAF_FILE_H
#define PANDORA_PAF_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"

class PafFile : public GenericFile {
private:
    void write_cluster(const Hits &cluster, uint32_t number_of_hits_in_read,
        const std::vector<MinimizerHitPtr> &sorted_minimizer_hits);

public:
    PafFile(const fs::path &filepath);
    void write_clusters(const MinimizerHitClusters &clusters);
};

#endif // PANDORA_PAF_FILE_H
