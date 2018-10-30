#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include <boost/log/trivial.hpp>

#include "minirecord.h"
#include "index.h"
#include "utils.h"


using namespace std;

Index::Index() = default;;

Index::~Index() {
    clear();
};

void Index::add_record(const uint64_t kmer, const uint32_t prg_id, const Path path, const uint32_t knode_id,
                       const bool strand) {
    //cout << "Add kmer " << kmer << " id, path, strand " << prg_id << ", " << path << ", " << strand << endl;
    auto it = minhash.find(kmer);
    if (it == minhash.end()) {
        auto *newv = new vector<MiniRecord>;
        newv->reserve(20);
        newv->emplace_back(MiniRecord(prg_id, path, knode_id, strand));
        minhash.insert(pair<uint64_t, vector<MiniRecord> *>(kmer, newv));
        //cout << "New minhash size: " << minhash.size() << endl; 
    } else {
        MiniRecord mr(prg_id, path, knode_id, strand);
        if (find(it->second->begin(), it->second->end(), mr) == it->second->end()) {
            it->second->push_back(mr);
        }
        //cout << "New minhash entry for  kmer " << kmer << endl;
    }
}

void Index::clear() {
    for (auto it = minhash.begin(); it != minhash.end();) {
        delete it->second;
        it = minhash.erase(it);
    }
}

void Index::save(const string &prgfile, uint32_t w, uint32_t k) {
    BOOST_LOG_TRIVIAL(debug) << "Saving index";
    ofstream handle;
    handle.open(prgfile + ".k" + to_string(k) + ".w" + to_string(w) + ".idx");

    handle << minhash.size() << endl;

    for (auto &it : minhash) {
        handle << it.first << "\t" << it.second->size();
        for (uint32_t j = 0; j != it.second->size(); ++j) {
            handle << "\t" << (*(it.second))[j];
        }
        handle << endl;

    }
    handle.close();
    BOOST_LOG_TRIVIAL(debug) << "Finished saving " << minhash.size() << " entries to file";
}

void Index::load(const string &prgfile, uint32_t w, uint32_t k) {
    BOOST_LOG_TRIVIAL(debug) << "Loading index";
    BOOST_LOG_TRIVIAL(debug) << "File is " << prgfile << ".k" << to_string(k) << ".w" << to_string(w) << ".idx";
    //string line;
    //vector<string> vstring;
    uint32_t key;
    size_t size;
    int c;
    MiniRecord mr;
    bool first = true;
    //vector<MiniRecord> vmr;

    ifstream myfile(prgfile + ".k" + to_string(k) + ".w" + to_string(w) + ".idx");
    if (myfile.is_open()) {
        while (myfile.good()) {
            c = myfile.peek();
            if (isdigit(c) and first) {
                myfile >> size;
                minhash.reserve(size);
                first = false;
                myfile.ignore(1, '\t');
            } else if (isdigit(c) and !first) {
                myfile >> key;
                myfile.ignore(1, '\t');
                myfile >> size;
                auto *vmr = new vector<MiniRecord>;
                vmr->reserve(size);
                minhash[key] = vmr;
                myfile.ignore(1, '\t');
            } else if (c == EOF) {
                break;
            } else {
                myfile >> mr;
                minhash[key]->push_back(mr);
                myfile.ignore(1, '\t');
            }
        }
    } else {
        cerr << "Unable to open index file " << prgfile << ".idx" << endl;
        exit(1);
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished loading " << minhash.size() << " entries to index";
}
	    

