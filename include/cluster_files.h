#ifndef PANDORA_CLUSTER_FILES_H
#define PANDORA_CLUSTER_FILES_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"

class ClusterDefFile : public GenericFile {
private:
public:
    ClusterDefFile(const fs::path &filepath, bool is_fake_file = false)
        : GenericFile(filepath, is_fake_file){
        (*this) << "read\tprg\tstatus\tcluster_size\tnb_of_repeated_mini\tnb_of_unique_mini\tlength_based_threshold\tmin_cluster_size\tread_hits_positions\tdistances_between_hits\n";
    }
};


class ClusterFilterFile : public GenericFile {
private:
public:
    ClusterFilterFile(const fs::path &filepath, bool is_fake_file = false)
        : GenericFile(filepath, is_fake_file){
        (*this) << "read\tprg\tnb_of_unique_minimisers\tstatus\t"
                   "removed_cluster_nb_of_unique_minimisers\tremoved_cluster_start\t"
                   "removed_cluster_end\tfavoured_cluster_start\t"
                   "favoured_cluster_end\toverlap\n";
    }
};

#endif // PANDORA_CLUSTER_FILES_H
