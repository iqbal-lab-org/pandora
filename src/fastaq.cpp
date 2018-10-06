#include <fstream>
#include <iostream>
#include <cassert>
#include <cctype>
#include <unordered_map>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "fastaq.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

Fastaq::Fastaq(bool gz, bool fq) : gzipped(gz), fastq(fq) {}

char Fastaq::covg_to_score(const uint_least16_t &covg, const uint_least16_t &global_covg) {
    if (2 * global_covg < covg) {
        cout << "Found a base with a coverage way too high, so giving it a score of 0" << endl;
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

void Fastaq::add_entry(const std::string &name,
                       const std::string &sequence,
                       const std::vector<uint32_t> &covgs,
                       const uint_least16_t global_covg,
                       const string header) {

    assert(name != "");
    assert(covgs.size() == sequence.length());
    assert(global_covg != 0);

    char score[covgs.size() + 1];
    auto i = 0;
    for (auto covg : covgs) {
        score[i] = covg_to_score(covg, global_covg);
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
                       const string header) {

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
    if (filepath.length() > 2 and filepath.substr(filepath.length() - 2) == "gz" and gzipped == false) {
        gzipped = true;
    } else if (filepath.length() > 2 and filepath.substr(filepath.length() - 2) != "gz" and gzipped == true) {
        gzipped = false;
    }
    ofstream file(filepath, ios_base::out | ios_base::binary | ios_base::trunc);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    if (gzipped) {
        out.push(boost::iostreams::gzip_compressor());
    }
    out.push(file);

    std::ostream outf(&out);
    outf << *this;
}


// returns coverage as just the number of reads in the fastaq
double Fastaq::calculate_coverage() const {
    return sequences.size();
}

// calculates coverage as number of bases / length of a given reference
double
Fastaq::calculate_kmer_coverage(const unsigned long &ref_length, const unsigned int k, const double &error_rate) const {
    const auto D{this->calculate_coverage()};

    return (D * (ref_length - k + 1)) / (ref_length * pow(1-error_rate, k));
}

bool Fastaq::operator==(const Fastaq &y) const {
    if (fastq != y.fastq) { return false; }
    if (names.size() != y.names.size()) { return false; }
    for (auto name : names) {
        if (find(y.names.begin(), y.names.end(), name) == y.names.end()) { return false; }
        if (y.sequences.find(name) == y.sequences.end()) { return false; }
        if (y.sequences.at(name) != sequences.at(name)) { return false; }
        if (fastq and y.scores.find(name) == y.scores.end()) { return false; }
        if (y.scores.at(name) != scores.at(name)) { return false; }
    }
    for (auto name : y.names) {
        if (find(names.begin(), names.end(), name) == names.end()) { return false; }
    }
    return true;
}

bool Fastaq::operator!=(const Fastaq &y) const {
    return !(*this == y);
}

std::ostream &operator<<(std::ostream &out, Fastaq const &data) {
    for (const auto name : data.names) {
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
    string name, seq, score, header;
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
