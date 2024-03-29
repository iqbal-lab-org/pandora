#ifndef PANDORA_PAF_FILE_H
#define PANDORA_PAF_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"
#include "seq.h"


class PafFile : public GenericFile {
private:
    const std::vector<std::string> &prg_names;
    void write_cluster(const Seq &seq, const MinimizerHits &cluster,
        const std::vector<MinimizerHitPtr> &all_sorted_minimizer_hits);

public:
    PafFile(const fs::path &filepath,
        const std::vector<std::string> &prg_names,
        bool is_fake_file = false);
    void write_clusters(const Seq &seq, const MinimizerHitClusters &clusters);
};

#endif // PANDORA_PAF_FILE_H
