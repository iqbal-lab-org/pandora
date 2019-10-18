#ifndef __FASTAQ_HANDLER_H_INCLUDED__ // if fastaq_handler.h hasn't been included yet...
#define __FASTAQ_HANDLER_H_INCLUDED__

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <cstdint>
#include <fstream>
#include <string>

namespace logging = boost::log;

struct FastaqHandler {
    bool gzipped;
    std::ifstream fastaq_file;
    boost::iostreams::filtering_istreambuf inbuf;
    std::istream instream;
    std::string line;
    std::string name;
    std::string read;
    uint32_t num_reads_parsed;

    FastaqHandler(const std::string&);

    ~FastaqHandler();

    bool eof();

    void get_next();

    void skip_next();

    void get_id(const uint32_t&);

    void close();
};

#endif
