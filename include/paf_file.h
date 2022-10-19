#ifndef PANDORA_PAF_FILE_H
#define PANDORA_PAF_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"
#include "seq.h"


class PafFile : public GenericFile {
private:
    const std::vector<std::shared_ptr<LocalPRG>>& prgs;
    void write_cluster(const Seq &seq, const Hits &cluster,
        const std::vector<MinimizerHitPtr> &all_sorted_minimizer_hits);

public:
    PafFile(const fs::path &filepath,
        const std::vector<std::shared_ptr<LocalPRG>>& prgs,
        bool is_fake_file = false);
    void write_clusters(const Seq &seq, const MinimizerHitClusters &clusters);
};

#endif // PANDORA_PAF_FILE_H
