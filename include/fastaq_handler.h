#ifndef __FASTAQ_HANDLER_H_INCLUDED__ // if fastaq_handler.h hasn't been included yet...
#define __FASTAQ_HANDLER_H_INCLUDED__

#include <string>
#include <cstdint>
#include <fstream>
#include <exception>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

namespace logging = boost::log;

struct FastaqHandler {
    const std::string &filepath;
    bool gzipped;
    gzFile fastaq_file;
    kseq_t* inbuf;
    std::string name;
    std::string read;
    uint32_t num_reads_parsed;
    int read_status;   // see https://github.com/attractivechaos/klib/blob/928581a78413bed4efa956731b35b18a638f20f3/kseq.h#L171
    int closed_status; // see https://github.com/attractivechaos/klib/blob/928581a78413bed4efa956731b35b18a638f20f3/kseq.h#L171


    FastaqHandler(const std::string&);

    ~FastaqHandler();

    bool eof() const;

    void get_next();

    void skip_next();

    void get_id(const uint32_t&);

    void close();

    bool is_closed() const;

};

#endif
