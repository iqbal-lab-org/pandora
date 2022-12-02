#ifndef PANDORA_GENERIC_FILE_H
#define PANDORA_GENERIC_FILE_H
#include <fstream>
#include <iostream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class GenericFile {
private:
    std::ofstream file_handler;
public:
    const bool is_fake_file;
    GenericFile(const fs::path &filepath, bool is_fake_file = false);

    virtual ~GenericFile() {
        if (not is_fake_file) {
            file_handler.close();
        }
    }

    template<class InputType>
    GenericFile& operator<<(const InputType &data) {
        if (not is_fake_file) {
            file_handler << data;
        }
        return *this;
    }
};

#endif
