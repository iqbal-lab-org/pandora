#ifndef __FASTAQ_HANDLER_H_INCLUDED__   // if fastaq_handler.h hasn't been included yet...
#define __FASTAQ_HANDLER_H_INCLUDED__

#include <string>
#include <cstdint>
//#include <iostream>
#include <fstream>

struct FastaqHandler {
    std::ifstream fastaq_file;
    std::string line;
    std::string name;
    std::string read;
    uint32_t num_reads_parsed;

    FastaqHandler(const std::string&);
    void get_next();
    void skip_next();
    void get_id(const uint32_t&);
    void close();
};

#endif
