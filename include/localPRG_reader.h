#ifndef PANDORA_LOCALPRG_READER_H
#define PANDORA_LOCALPRG_READER_H

#include "localPRG.h"
#include "fastaq_handler.h"
#include "forward_declarations.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class LocalPRGReaderGenerator
{
public:
    typedef std::shared_ptr<LocalPRG> result_type;
    LocalPRGReaderGenerator(const fs::path &filepath) : fh(filepath.string()), id(0) { }
    std::shared_ptr<LocalPRG> operator()();
private:
    FastaqHandler fh;
    uint32_t id;
};

class LocalPRGReader {
public:
    static std::vector<std::shared_ptr<LocalPRG>> read_prg_file_as_vector(const fs::path &filepath);

    /** Read PRGs and return generator
     * @return We need to return both the generator and iterator because the lifetime
     * of the generator should be the same as the iterator. Otherwise the FastaqHandler
     * of the generator will be closed as soon as we return.
     */
    static inline LocalPRGGeneratorAndIterator read_prg_file_as_generator(const fs::path &filepath) {
        LocalPRGReaderGeneratorPtr local_prg_reader = std::make_shared<LocalPRGReaderGenerator>(filepath);
        LocalPRGReaderGeneratorIterator prg_it(boost::make_generator_iterator(*local_prg_reader));
        return std::make_pair(local_prg_reader, prg_it);
    }
};

#endif // PANDORA_LOCALPRG_READER_H
