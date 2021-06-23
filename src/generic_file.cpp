#include "generic_file.h"
#include <iostream>
#include "utils.h"

GenericFile::GenericFile(const fs::path &filepath) {
    open_file_for_writing(filepath.string(), file_handler);
}
GenericFile::~GenericFile() {
    file_handler.close();
}