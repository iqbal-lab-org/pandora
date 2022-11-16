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

    std::string get_file_given_path(const std::string &zip_path) {
        zip_int64_t zip_file_location = zip_name_locate(archive, zip_path.c_str(), 0);
        zip_file *pFile = zip_fopen_index(archive, zip_file_location, 0);
        /*
        if( NULL == pFile){
            qDebug()<<"fopen is NULL"<<zip_strerror(m_zip);
            return QByteArray();
        }
         */
        struct zip_stat stat;
        int rStat = zip_stat_index(archive, zip_file_location, 0, &stat);
        /*
        if( -1 == rStat ){
            qDebug()<<"stat failed : "<<zip_strerror(m_zip);
            return QByteArray();
        }
         */
        const int length = stat.size;
        char buffer[length +1 ];
        int rRead = zip_fread(pFile,buffer,sizeof(buffer));
        /*
        if( -1 == rRead ){
            qDebug()<<"read failed : "<<zip_strerror(m_zip);
            return QByteArray();
        }*/
        return std::string(buffer);
    }
};

#endif // PANDORA_ZIP_FILE_H
