#include "zip_file.h"

std::pair<zip_file*, struct zip_stat> ZipFileReader::open_file_inside_zip(
    struct zip *archive, const std::string &zip_path) {
    zip_int64_t zip_file_location = zip_name_locate(archive, zip_path.c_str(), 0);
    zip_file* zipsub_file;
    if ((zipsub_file =zip_fopen_index(archive, zip_file_location, 0)) == nullptr) {
        fatal_error("Unable to open file ", zip_path, " inside zip file for reading");
    }

    struct zip_stat stat;
    if (zip_stat_index(archive, zip_file_location, 0, &stat) == -1) {
        fatal_error("Unable to open file ", zip_path, " inside zip file for reading");
    }

    return std::make_pair(zipsub_file, stat);
}


std::vector<std::string> ZipFileReader::read_full_text_file(const std::string &zip_path) {
    zip_file* zipsub_file;
    struct zip_stat stat;
    std::tie(zipsub_file, stat) = open_file_inside_zip(archive, zip_path);

    const int length = stat.size;
    char buffer[length +1 ];
    if(zip_fread(zipsub_file, buffer, length) != length){
        fatal_error("Unable to read from: ", zip_path);
    }
    buffer[length] = '\0';

    if(zip_fclose(zipsub_file) != 0) {
        fatal_error("Unable to close: ", zip_path);
    }

    return split(buffer, "\n");
}
