#include <fstream>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "fastaq.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

Fastaq::Fastaq(bool gz, bool fq) : gzipped(gz), fastq(fq) {}

char Fastaq::covg_to_score(const uint_least16_t& covg, const uint_least16_t& global_covg){
    assert(global_covg >= covg);
    int c = 40*covg/global_covg + 33;
    char ascii_c = static_cast<char>(c);
    return ascii_c;
}

void Fastaq::add_entry(const std::string & name,
                  const std::string & sequence,
                  const std::vector <uint32_t> & covgs,
                  const uint_least16_t global_covg){

    assert(name != "");
    assert(covgs.size() == sequence.length());
    assert(global_covg!=0);

    char score[covgs.size()];
    auto i = 0;
    for (auto covg : covgs){
        score[i] = covg_to_score(covg,global_covg);
        i++;
    }

    names.push_back(name);
    sequences[name] = sequence;
    scores[name] = score;
}

void Fastaq::save(const std::string & filepath) {
    ofstream file(filepath, ios_base::out | ios_base::binary | ios_base::trunc);
    boost::iostreams::filtering_streambuf<boost::iostreams::output> out;
    if (gzipped)
        out.push(boost::iostreams::gzip_compressor());
    out.push(file);

    std::ostream outf(&out);
    outf << *this;
}

bool Fastaq::operator==(const Fastaq &y) const {
    if (fastq != y.fastq) { return false;}
    if (names.size() != y.names.size()) { return false;}
    for (auto name : names) {
        if (find(y.names.begin(), y.names.end(), name) == y.names.end()) {return false;}
        if (y.sequences.find(name) == y.sequences.end()) { return false;}
        if (y.sequences.at(name) != sequences.at(name)) {return false;}
        if (fastq and y.scores.find(name) == y.scores.end()) { return false;}
        if (y.scores.at(name) != scores.at(name)) {return false;}
    }
    for (auto name : y.names) {
        if (find(names.begin(), names.end(), name) == names.end()) { return false; }
    }
    return true;
}

bool Fastaq::operator!=(const Fastaq &y) const {
    return !(*this==y);
}

std::ostream &operator<<(std::ostream &out, Fastaq const &data) {
    for (const auto name : data.names){
        if (data.fastq)
            out << "@";
        else
            out << ">";
        out << name << "\n";
        out << data.sequences.at(name) << "\n";
        if (data.fastq) {
            out << "+\n";
            out << data.scores.at(name) << "\n";
        }
    }
    return out;
}

std::istream &operator>>(std::istream &in, Fastaq &data) {
    string name, seq, score;
    /*in >> name;
    in >> seq;
    data.names.push_back(name);
    data.sequences[name] = seq;*/

    char c = in.peek();
    while (c != EOF) {
        in.ignore(1, '>');
        in >> name;
        in >> seq;
        data.names.push_back(name);
        data.sequences[name] = seq;
        in.ignore(1, '\n');
        c = in.peek();
        cout << c << " " << (c == '+') << endl;
        if (c == '+') {
            data.fastq = true;
            in.ignore(1, '+');
            in >> score;
            data.scores[name] = score;
        }
        c = in.peek();
    }

    return in;
}