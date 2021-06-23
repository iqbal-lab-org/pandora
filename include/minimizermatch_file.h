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
    const std::vector<std::shared_ptr<LocalPRG>>& prgs;
public:
    MinimizerMatchFile(const fs::path &filepath,
                       // just to convert prg IDs to prg names
                       const std::vector<std::shared_ptr<LocalPRG>>& prgs);
    void write_hits(const Seq &seq, const Hits &hits);
};


#endif // PANDORA_MINMATCH_FILE_H
