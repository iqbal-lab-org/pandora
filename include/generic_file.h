#ifndef PANDORA_GENERIC_FILE_H
#define PANDORA_GENERIC_FILE_H
#include <fstream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class GenericFile {
protected:
    std::ofstream file_handler;
public:
    GenericFile(const fs::path &filepath);
    virtual ~GenericFile();
};

#endif
