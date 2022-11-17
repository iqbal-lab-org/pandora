#ifndef PANDORA_ZIP_FILE_H
#define PANDORA_ZIP_FILE_H

#include <archive.h>
#include <archive_entry.h>
#include <zip.h>
#include "utils.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class ZipFileWriter {
private:
    struct archive *zip_archive;
    struct archive_entry *current_zip_archive_entry;

    inline void cleanup_previous_entry() {
        if (current_zip_archive_entry) {
            archive_entry_free(current_zip_archive_entry);
            current_zip_archive_entry = nullptr;
        }
    }

public:
    // note: not explicit as we want to build from everything we can convert to a path
    ZipFileWriter(const fs::path &path) : current_zip_archive_entry(nullptr) {
        zip_archive = archive_write_new();
        archive_write_set_format_zip(zip_archive);
        archive_write_open_filename(zip_archive, path.string().c_str());
    }
    virtual ~ZipFileWriter() {
        cleanup_previous_entry();
        archive_write_close(zip_archive);
        archive_write_free(zip_archive);
    }
    ZipFileWriter(ZipFileWriter&& other) = default; // move default constructor
    ZipFileWriter& operator=(ZipFileWriter&& other) = delete; // move assignment operator

    // not allowed to copy/assign
    ZipFileWriter(const ZipFileWriter& other) = delete;
    ZipFileWriter& operator=(const ZipFileWriter& other) = delete;


    inline void prepare_new_entry(const std::string &zip_path) {
        cleanup_previous_entry();
        current_zip_archive_entry = archive_entry_new();
        archive_entry_set_mode(current_zip_archive_entry, AE_IFREG | 0644);
        archive_entry_set_pathname(current_zip_archive_entry, zip_path.c_str());
        if (archive_write_header(zip_archive, current_zip_archive_entry) != ARCHIVE_OK) {
            fatal_error("Error writing zip archive entry");
        }
    }

    inline void write_data(const std::string &data) {
        if (archive_write_data(zip_archive, data.c_str(), data.size()) != data.size()) {
            fatal_error("Error writing zip data");
        }
    }
};

class ZipFileReader {
private:
    struct zip *archive;
    std::string read_full_text_file_as_single_string(const std::string &zip_path);
    inline std::vector<std::string> read_full_text_file(const std::string &zip_path) {
        return split(read_full_text_file_as_single_string(zip_path), "\n");
    }

    static std::pair<zip_file*, struct zip_stat> open_file_inside_zip(
        struct zip *archive, const std::string &zip_path);
public:
    ZipFileReader(const fs::path &path) {
        int err;
        if ((archive = zip_open(path.c_str(), ZIP_RDONLY, &err)) == nullptr) {
            fatal_error("Cannot open zip file ", path.string());
        }
    }

    ~ZipFileReader() {
        if (zip_close(archive) == -1) {
            fatal_error("Can't close zip archive");
        }
    }

    inline std::pair<uint32_t, uint32_t> read_w_and_k() {
        std::vector<std::string> w_and_k_as_strings = read_full_text_file("_metadata");
        return std::make_pair(std::stoi(w_and_k_as_strings[0]),
                              std::stoi(w_and_k_as_strings[1]));
    }

    inline std::vector<std::string> read_prg_names() {
        return read_full_text_file("_prg_names");
    }

    inline std::vector<uint32_t> read_prg_min_path_lengths() {
        std::vector<std::string> prg_min_path_lengths_as_strs = read_full_text_file("_prg_min_path_lengths");
        std::vector<uint32_t> prg_min_path_lengths(prg_min_path_lengths_as_strs.size());
        std::transform(prg_min_path_lengths_as_strs.begin(),
                       prg_min_path_lengths_as_strs.end(),
                       prg_min_path_lengths.begin(),
                       [](const std::string &str) { return std::stoi(str); }
        );
        return prg_min_path_lengths;
    }

    friend class ZipInstreamBuffer;
    friend class Index;
};


// Adapted from https://stackoverflow.com/a/9140328/5264075
class ZipInstreamBuffer : public std::streambuf {
private:
    zip_file* file_;
    enum { s_size = 8196 };
    char buffer_[s_size];

public:
    ZipInstreamBuffer(struct zip *archive, const std::string &zip_path){
        struct zip_stat stat;
        std::tie(file_, stat) = ZipFileReader::open_file_inside_zip(archive, zip_path);
    }
    ~ZipInstreamBuffer(){
        if(zip_fclose(file_) != 0) {
            fatal_error("Unable to close");
        }
    }

protected:
    int underflow() override {
        int rc(zip_fread(this->file_, this->buffer_, s_size));
        this->setg(this->buffer_, this->buffer_,
            this->buffer_ + std::max(0, rc));
        return this->gptr() == this->egptr()
            ? traits_type::eof()
            : traits_type::to_int_type(*this->gptr());
    }
};

#endif // PANDORA_ZIP_FILE_H
