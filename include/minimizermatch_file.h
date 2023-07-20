#ifndef PANDORA_MINMATCH_FILE_H
#define PANDORA_MINMATCH_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"
#include "seq.h"
#include "localPRG.h"

class MinimizerMatchFile : public GenericFile {
private:
    const std::vector<std::string> &prg_names;
public:
    MinimizerMatchFile(const fs::path &filepath, const std::vector<std::string> &prg_names,
                       bool is_fake_file = false);
    void write_hits(const Seq &seq, const MinimizerHits &hits, const uint32_t k);
};


#endif // PANDORA_MINMATCH_FILE_H
