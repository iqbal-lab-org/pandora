#include <string>
#include <iostream>
#include "fastaq_handler.h"

FastaqHandler::FastaqHandler(const std::string filepath)
    : closed(false)
    , filepath(filepath)
    , num_reads_parsed(0)
{
    this->fastaq_file = gzopen(filepath.c_str(), "r");
    if (this->fastaq_file == nullptr) {
        throw std::ios_base::failure("Unable to open " + this->filepath);
    }
    this->inbuf = kseq_init(this->fastaq_file);
}

FastaqHandler::~FastaqHandler() { this->close(); }

bool FastaqHandler::eof() const { return ks_eof(this->inbuf->f); }

void FastaqHandler::get_next()
{
    if (this->eof()) {
        throw std::out_of_range("Read requested after the end of file was reached");
    }
    int read_status = kseq_read(this->inbuf);

    // if not eof but we get -1 here then it was an empty file/read/line
    if (read_status == -1) {
        throw std::out_of_range("Read requested after the end of file was reached");
    }
    if (read_status == -2) {
        throw std::runtime_error("Truncated quality string detected");
    } else if (read_status == -3) {
        throw std::ios_base::failure("Error reading " + this->filepath);
    }

    ++this->num_reads_parsed;
    this->name = this->inbuf->name.s;
    this->read = this->inbuf->seq.s;
}

void FastaqHandler::get_nth_read(const uint32_t& idx)
{
    // edge case where no reads have been loaded yet
    if (this->num_reads_parsed == 0) {
        this->get_next();
    }
    const uint32_t one_based_idx = idx + 1;
    if (one_based_idx < this->num_reads_parsed) {
        num_reads_parsed = 0;
        name.clear();
        read.clear();
        gzrewind(this->fastaq_file);
        kseq_rewind(this->inbuf);
    }

    while (this->num_reads_parsed < one_based_idx) {
        try {
            this->get_next();
        } catch (std::out_of_range& err) {
            throw std::out_of_range("Requested a read past the end of the file.");
        }
    }
}

void FastaqHandler::close()
{
    if (!this->is_closed()) {
        const auto closed_status = gzclose(this->fastaq_file);
        kseq_destroy(this->inbuf);

        if (closed_status != Z_OK) {
            std::ostringstream err_msg;
            err_msg << "Failed to close " << this->filepath
                    << ". Got zlib return code: " << closed_status << std::endl;
            throw std::ios_base::failure(err_msg.str());
        }
        this->closed = true;
    }
}

bool FastaqHandler::is_closed() const { return this->closed; }
