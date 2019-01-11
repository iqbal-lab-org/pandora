#include <fstream>
#include <iostream>
#include <cassert>
#include <cctype>
#include <unordered_map>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/log/trivial.hpp>

#include "fastaq.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


Fastaq::Fastaq(bool gz, bool fq) : gzipped(gz), fastq(fq) {}

char Fastaq::covg_to_score(const uint_least16_t &covg, const uint_least16_t &global_covg, const bool &alt) {
    // this alternate conversion maps coverage ASCII of value up to 93 and then anything over this is just mapped to
    // the ASCII value for 93
    if (alt) {
        return Fastaq::alt_covg_to_score(covg);
    }
    // Rachel's original (and default) coverage to ASCII conversion function
    if (2 * global_covg < covg) {
        BOOST_LOG_TRIVIAL(warning) << "Found a base with a coverage way too high, so giving it a score of 0";
        return '!';
    }

    int c;
    if (global_covg >= covg) {
        c = 40 * covg / global_covg + 33;
    } else {
        c = 40 * (2 * global_covg - covg) / global_covg + 33;
    }
    char ascii_c = static_cast<char>(c);
    return ascii_c;
}


char Fastaq::alt_covg_to_score(const uint_least16_t &covg) {
    // use ASCII chars 33 - 126 as these are the printable ones
    const uint_least16_t max{126 - 33};
    uint_least16_t ascii_val;

    if (covg > max) {  // coverage is outside the range of printable ASCIIs
        ascii_val = 126;
    } else {
        ascii_val = covg + 33;
    }
    return static_cast<char>(ascii_val);
}


void Fastaq::add_entry(const std::string &name,
                       const std::string &sequence,
                       const std::vector<uint32_t> &covgs,
                       const uint_least16_t global_covg,
                       const std::string header) {

    assert(name != "");
    assert(covgs.size() == sequence.length());
    auto mod_global_covg = global_covg;
    if (global_covg < 1) {
        mod_global_covg = 1;
    }

    char score[covgs.size() + 1];
    auto i = 0;
    const bool alt_covg_conversion{false};
    for (const auto &covg: covgs) {
        score[i] = covg_to_score(covg, mod_global_covg, alt_covg_conversion);
        i++;
    }
    score[covgs.size()] = '\0';

    names.push_back(name);
    headers[name] = header;
    sequences[name] = sequence;
    scores[name] = score;
}

void Fastaq::add_entry(const std::string &name,
                       const std::string &sequence,
                       const std::string header) {

    assert(name != "");

    names.push_back(name);
    headers[name] = header;
    sequences[name] = sequence;
    scores[name] = {};
}

void Fastaq::clear() {
    names.clear();
    headers.clear();
    sequences.clear();
    scores.clear();
}

void Fastaq::save(const std::string &filepath) {
    if (filepath.length() > 2 and filepath.substr(filepath.length() - 2) == "gz" and !gzipped) {
        gzipped = true;
    } else if (filepath.length() > 2 and filepath.substr(filepath.length() - 2) != "gz" and gzipped) {
        gzipped = false;
    }
    std::ofstream file(filepath, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    if (gzipped) {
        out.push(boost::iostreams::gzip_compressor());
    }
    out.push(file);

    std::ostream outf(&out);
    outf << *this;
}


bool Fastaq::operator==(const Fastaq &y) const {
    if (fastq != y.fastq) { return false; }
    if (names.size() != y.names.size()) { return false; }
    for (const auto &name: names) {
        if (find(y.names.begin(), y.names.end(), name) == y.names.end()) { return false; }
        if (y.sequences.find(name) == y.sequences.end()) { return false; }
        if (y.sequences.at(name) != sequences.at(name)) { return false; }
        if (fastq and y.scores.find(name) == y.scores.end()) { return false; }
        if (y.scores.at(name) != scores.at(name)) { return false; }
    }
    for (const auto &name: y.names) {
        if (find(names.begin(), names.end(), name) == names.end()) { return false; }
    }
    return true;
}

bool Fastaq::operator!=(const Fastaq &y) const {
    return !(*this == y);
}

std::ostream &operator<<(std::ostream &out, Fastaq const &data) {
    for (const auto &name : data.names) {
        if (data.fastq) {
            out << "@";
        } else {
            out << ">";
        }
        out << name;
        if (data.headers.at(name) != "") {
            out << data.headers.at(name);
        }
        out << "\n";
        out << data.sequences.at(name) << "\n";
        if (data.fastq) {
            out << "+\n";
            out << data.scores.at(name) << "\n";
        }
    }
    return out;
}

std::istream &operator>>(std::istream &in, Fastaq &data) {
    std::string name, seq, score, header;
    /*in >> name;
    in >> seq;
    data.names.push_back(name);
    data.sequences[name] = seq;*/

    char c = in.peek();
    while (c != EOF) {

        in.ignore(1, '>');
        in >> name;
        data.names.push_back(name);

        data.headers[name] = "";
        c = in.peek();
        while (isspace(c) and c != '\n') {
            in >> header;
            data.headers[name] += " " + header;
            c = in.peek();
        }

        in >> seq;
        data.sequences[name] = seq;

        in.ignore(1, '\n');
        c = in.peek();
        if (c == '+') {
            data.fastq = true;
            in.ignore(1, '+');
            in >> score;
            data.scores[name] = score;
            c = in.peek();
        }
    }

    return in;
}
