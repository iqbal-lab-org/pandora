#include "zip_file.h"

void ZipFile::prepare_new_entry(const std::string &zip_path) {
    cleanup_previous_entry();
    current_zip_archive_entry = archive_entry_new();
    archive_entry_set_mode(current_zip_archive_entry, AE_IFREG | 0644);
    archive_entry_set_pathname(current_zip_archive_entry, zip_path.c_str());
    if (archive_write_header(zip_archive, current_zip_archive_entry) != ARCHIVE_OK) {
        fatal_error("Error writing zip archive entry");
    }
}
