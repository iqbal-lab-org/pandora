#include <string>
#include <cstring>
#include <iostream>
//#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "fastaq_handler.h"
#include "utils.h"

FastaqHandler::FastaqHandler(const std::string& filepath)
    : gzipped(false)
    , instream(&inbuf)
    , num_reads_parsed(0)
{
    // level for boost logging
    //    logging::core::get()->set_filter(logging::trivial::severity >= g_log_level);

    BOOST_LOG_TRIVIAL(debug) << "Open fastaq file" << filepath;
    fastaq_file.open(filepath);
    if (not fastaq_file.is_open()) {
        std::cerr << "Unable to open fastaq file " << filepath << std::endl;
        std::exit(EXIT_FAILURE);
    }
    try {
        if (filepath.substr(filepath.length() - 2) == "gz") {
            inbuf.push(boost::iostreams::gzip_decompressor());
            gzipped = true;
        }
        inbuf.push(fastaq_file);
    } catch (const boost::iostreams::gzip_error& e) {
        std::cerr << "Problem transfering file contents to boost stream: " << e.what()
                  << '\n';
    }
}

FastaqHandler::~FastaqHandler() { close(); }

bool FastaqHandler::eof()
{
    int c = instream.peek();
    return (c == EOF);
}

void FastaqHandler::get_next()
{
    // cout << "next ";
    if (!line.empty() and (line[0] == '>' or line[0] == '@')) {
        // cout << "read name line " << num_reads_parsed << " " << line << endl;
        name = line.substr(1);
        ++num_reads_parsed;
        read.clear();
    }

    while (getline(instream, line).good()) {
        if (!line.empty() and line[0] == '+') {
            // skip this line and the qual score line
            getline(instream, line);
            // cout << "qual score line ." << line << "." << endl;
        } else if (line.empty() || line[0] == '>' || line[0] == '@') {
            if (!read.empty()
                or line.empty()) // ok we'll allow reads with no name, removed
            {
                return;
            }
            // cout << num_reads_parsed << " " << line << endl;
            name = line.substr(1);
            ++num_reads_parsed;
            read.clear();
        } else {
            // cout << "read line ." << line << "." << endl;
            read += line;
        }
    }
}

void FastaqHandler::skip_next()
{
    // cout << "skip ";
    if (!line.empty() and (line[0] == '>' or line[0] == '@')) {
        ++num_reads_parsed;
    }

    while (getline(instream, line).good()) {
        if (!line.empty() and line[0] == '+') {
            // skip this line and the qual score line
            getline(instream, line);
            // cout << "qual score line ." << line << "." << endl;
        } else if (line.empty() or line[0] == '>' or line[0] == '@') {
            return;
        }
    }
}

void print(std::ifstream& infile)
{
    char file;
    std::vector<char> read;
    uint i = 0;

    // Read infile to vector
    while (!infile.eof()) {
        infile >> file;
        read.push_back(file);
    }

    // Print read vector
    for (i = 0; i < read.size(); i++) {
        std::cout << read[i];
    }
}

void print(std::istream& infile)
{
    char file;
    std::vector<char> read;
    int i = 0;

    // Read infile to vector
    while (!infile.eof()) {
        infile >> file;
        read.push_back(file);
    }

    // Print read vector
    for (auto i = 0; i < read.size(); i++) {
        std::cout << read[i];
    }
}

void FastaqHandler::get_id(const uint32_t& id)
{
    const uint32_t one_based_id = id + 1;
    if (one_based_id < num_reads_parsed) {
        BOOST_LOG_TRIVIAL(warning)
            << "restart buffer as have id " << num_reads_parsed << " and want id "
            << one_based_id << " (" << id << ") with 0-based indexing.";
        num_reads_parsed = 0;
        name.clear();
        read.clear();
        line.clear();
        assert(name.empty());
        assert(read.empty());
        assert(line.empty());

        instream.ignore(std::numeric_limits<std::streamsize>::max());
        fastaq_file.seekg(0, fastaq_file.beg);
        instream.clear();
        inbuf.pop();
        if (gzipped) {
            inbuf.pop();
            inbuf.push(boost::iostreams::gzip_decompressor());
        }
        inbuf.push(fastaq_file);
    }

    while (id > 1 and num_reads_parsed < id) {
        skip_next();
        if (eof()) {
            break;
        }
    }

    while (num_reads_parsed <= id) {
        get_next();
        if (eof()) {
            break;
        }
    }
}

void FastaqHandler::close()
{
    BOOST_LOG_TRIVIAL(debug) << "Close fastaq file";
    fastaq_file.close();
}
