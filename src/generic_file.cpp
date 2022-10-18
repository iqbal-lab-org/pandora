#include "generic_file.h"
#include "utils.h"

GenericFile::GenericFile(const fs::path &filepath, bool is_fake_file)
  : is_fake_file(is_fake_file) {
    if (not is_fake_file) {
        open_file_for_writing(filepath.string(), file_handler);
    }
}
