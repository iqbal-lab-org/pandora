#ifndef PANDORA_ZIP_FILE_H
#define PANDORA_ZIP_FILE_H

#include <archive.h>
#include <archive_entry.h>
#include <zip.h>
#include "utils.h"
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;


constexpr size_t zip_file_buffer_size = 16*1024; //16 kb


// Use this class to write to Zip files
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
    // Note: this constructor is not explicit as we want to build from everything we can convert to a path
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

    /**
     * Creates a new entry in the zip file.
     * Data to this entry can be written through the write_data() method.
     * @param zip_path : Path of the entry in the zip file
     */
    inline void prepare_new_entry(const std::string &zip_path) {
        cleanup_previous_entry();
        current_zip_archive_entry = archive_entry_new();
        archive_entry_set_mode(current_zip_archive_entry, AE_IFREG | 0644);
        archive_entry_set_pathname(current_zip_archive_entry, zip_path.c_str());
        if (archive_write_header(zip_archive, current_zip_archive_entry) != ARCHIVE_OK) {
            fatal_error("Error writing zip archive entry");
        }
    }

    /**
     * Writes the given data to the prepared new entry done with prepare_new_entry() method.
     * Be careful when using this method, as it requires the whole data to be in memory
     * before writing. Can be used to write data in chunks.
     */
    void write_data(const std::string &data);

    /**
     * Writes the lines from the given ifstream to the prepared new entry done with prepare_new_entry() method.
     * Note that this just works for text files.
     * Use this to write text files from disk to this zip file.
     */
    inline void write_from_text_stream(std::ifstream& data_source) {
        std::string line;
        while (std::getline(data_source, line)) {
            line += "\n";
            write_data(line);
        }
    }
};


// Use this class to read from Zip files.
// There are some convenience methods here.
// See ZipIfstream class for reading from files from Zip using streams.
class ZipFileReader {
private:
    struct zip *archive;

    std::string read_full_text_file_as_single_string(const std::string &zip_path);
    inline std::vector<std::string> read_full_text_file(const std::string &zip_path) {
        return split(read_full_text_file_as_single_string(zip_path), "\n");
    }

    static std::pair<zip_file*, struct zip_stat> open_file_inside_zip(
        struct zip *archive, const std::string &zip_path);

    /**
     * Extract a file from the zip archive into a temp file
     * @param zip_path : The file to be extracted
     * @return A path to the temp file extracted. The caller is responsible for deleting this file.
     */
    fs::path extract_text_file(const std::string &zip_path);

public:
    // Note: this constructor is not explicit as we want to build from everything we can convert to a path
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

    ////////////////////////////////////////////////////////////////////////////////////
    // List of convenient read methods for pandora
    inline std::pair<uint32_t, uint32_t> read_w_and_k() {
        std::vector<std::string> w_and_k_as_strings = read_full_text_file("_metadata");
        return std::make_pair(std::stoi(w_and_k_as_strings[0]),
                              std::stoi(w_and_k_as_strings[1]));
    }
    inline std::vector<std::string> read_prg_names() {
        return read_full_text_file("_prg_names");
    }
    inline std::vector<uint32_t> read_int_values(const std::string &zip_path) {
        std::vector<std::string> values_as_strs = read_full_text_file(zip_path);
        std::vector<uint32_t> values(values_as_strs.size());
        std::transform(values_as_strs.begin(),
            values_as_strs.end(),
            values.begin(),
            [](const std::string &str) { return std::stoi(str); }
        );
        return values;
    }
    inline std::vector<uint32_t> read_prg_lengths() {
        return read_int_values("_prg_lengths");
    }
    inline std::vector<uint32_t> read_prg_min_path_lengths() {
        return read_int_values("_prg_min_path_lengths");
    }
    inline fs::path extract_prgs() {
        return extract_text_file("_prgs");
    }
    ////////////////////////////////////////////////////////////////////////////////////

    friend class ZipIfstream;
    friend class Index;
};

class ZipIfstream {
private:
    // Adapted from https://stackoverflow.com/a/9140328/5264075
    class ZipInstreamBuffer : public std::streambuf {
    private:
        zip_file* file_;
        char buffer_[zip_file_buffer_size];

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
            int rc(zip_fread(this->file_, this->buffer_, zip_file_buffer_size));
            this->setg(this->buffer_, this->buffer_,
                this->buffer_ + std::max(0, rc));
            return this->gptr() == this->egptr()
                ? traits_type::eof()
                : traits_type::to_int_type(*this->gptr());
        }
    };

    ZipInstreamBuffer zip_buffer;
    std::istream zip_ifs;
public:
    ZipIfstream(struct zip *archive, const std::string &zip_path) :
        zip_buffer(archive, zip_path), zip_ifs(&zip_buffer) { }

    inline auto good() const {
        return zip_ifs.good();
    }

    inline auto peek() {
        return zip_ifs.peek();
    }

    inline void ignore(std::streamsize size, char delim) {
        zip_ifs.ignore(size, delim);
    }

    template <class T>
    ZipIfstream& operator>>(T& data) {
        zip_ifs >> data;
        return *this;
    }

    std::string getline() {
        std::string line;
        std::getline(zip_ifs, line);
        return line;
    }
};




#endif // PANDORA_ZIP_FILE_H
