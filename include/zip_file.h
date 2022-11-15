#ifndef PANDORA_ZIP_FILE_H
#define PANDORA_ZIP_FILE_H

#include <archive.h>
#include <archive_entry.h>
#include "utils.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class ZipFile {
private:
    struct archive *zip_archive;
    struct archive_entry *current_zip_archive_entry;

    inline void cleanup_previous_entry() {
        if (current_zip_archive_entry) {
            archive_entry_free(current_zip_archive_entry);
            current_zip_archive_entry = nullptr;
        }
    }
    inline void write_data_core(void* data, size_t data_size) {
        if (archive_write_data(zip_archive, data, data_size) != data_size) {
            fatal_error("Error writing zip data");
        }
    }

public:
    // note: not explicit as we want to build from everything we can convert to a path
    ZipFile(const fs::path &path) : current_zip_archive_entry(nullptr) {
        zip_archive = archive_write_new();
        archive_write_set_format_zip(zip_archive);
        archive_write_open_filename(zip_archive, path.string().c_str());
    }
    virtual ~ZipFile() {
        cleanup_previous_entry();
        archive_write_close(zip_archive);
        archive_write_free(zip_archive);
    }
    ZipFile(ZipFile&& other) = default; // move default constructor
    ZipFile& operator=(ZipFile&& other) = delete; // move assignment operator

    // not allowed to copy/assign
    ZipFile(const ZipFile& other) = delete;
    ZipFile& operator=(const ZipFile& other) = delete;


    void prepare_new_entry(const std::string &zip_path);
    inline void write_data(const std::string &data, bool synchronize) {
        if (synchronize) {
#pragma omp critical(pandora_zip_file__write_data)
            {
                write_data_core((void*)data.c_str(), data.size());
            }
        } else {
            write_data_core((void*)data.c_str(), data.size());
        }
    }
};

#endif // PANDORA_ZIP_FILE_H
