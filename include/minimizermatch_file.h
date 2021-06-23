#ifndef PANDORA_MINMATCH_FILE_H
#define PANDORA_MINMATCH_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"

class MinimizerMatchFile : public GenericFile {
public:
    MinimizerMatchFile(const fs::path &filepath);
    void write_hits(const Hits &hits);
};


#endif // PANDORA_MINMATCH_FILE_H
