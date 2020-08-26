#include <string>
#include <iostream>
#include "fastaq_handler.h"

FastaqHandler::FastaqHandler(const std::string filepath)
    : closed_status(3)
    , filepath(filepath)
    , num_reads_parsed(0)
{
    this->fastaq_file = gzopen(filepath.c_str(), "r");
    if (this->fastaq_file == nullptr) {
        throw "Unable to open " + this->filepath;
    }
    this->inbuf = kseq_init(this->fastaq_file);
}

FastaqHandler::~FastaqHandler() { this->close(); }

bool FastaqHandler::eof() const { return ks_eof(this->inbuf->f); }

void FastaqHandler::get_next()
{
    int read_status = kseq_read(this->inbuf);

    bool no_more_reads_available = read_status == -1;
    if (no_more_reads_available) {
        return;
    }

    if (read_status == -2) {
        throw "Truncated quality string detected";
    } else if (read_status == -3) {
        throw "Error reading " + this->filepath;
    }

    ++this->num_reads_parsed;
    this->name = this->inbuf->name.s;
    this->read = this->inbuf->seq.s;
}

void FastaqHandler::get_id(const uint32_t& id)
{
    const uint32_t one_based_id = id + 1;
    if (one_based_id < this->num_reads_parsed) {
        num_reads_parsed = 0;
        name.clear();
        read.clear();
        gzrewind(this->fastaq_file);
        kseq_rewind(this->inbuf);
    }

    while (id > 1 and num_reads_parsed < id) {
        get_next();
        if (eof()) {
            break;
        }
    }

    while (num_reads_parsed <= id) {
        get_next();
        if (eof() && num_reads_parsed < id) {
            throw std::out_of_range("Requested a read past the end of the file.");
        }
    }
}

void FastaqHandler::close()
{
    if (!this->is_closed()) {
        this->closed_status = gzclose(this->fastaq_file);
    }
}

bool FastaqHandler::is_closed() const { return this->closed_status == Z_OK; }
