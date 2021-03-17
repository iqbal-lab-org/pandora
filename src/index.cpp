#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <algorithm>

#include <boost/log/trivial.hpp>

#include "minirecord.h"
#include "index.h"
#include "localPRG.h"

/**
 * Adds a k-mer to the index. This is *just* called to add minimizers.
 *
 * @param kmer : the hash value of the minimizer (smallest hash of the canonicals)
 * @param prg_id : the prg from where this k-mer comes from
 * @param path : the path of this kmer in the prg
 * @param knode_id : the id of the node representing this kmer in the KmerGraph
 * @param strand : the strand
 */
void Index::add_record(const uint64_t kmer, const uint32_t prg_id,
    const prg::Path& path, const uint32_t knode_id, const bool strand)
{
    auto it = minhash.find(kmer); // checks if kmer is in minhash
    if (it == minhash.end()) { // no
        auto* newv = new std::vector<MiniRecord>; // get a new vector of MiniRecords -
                                                  // TODO: is this a memory leak?
        newv->reserve(20);
        newv->emplace_back(MiniRecord(prg_id, path, knode_id, strand));
        minhash.insert(std::pair<uint64_t, std::vector<MiniRecord>*>(
            kmer, newv)); // add this record to minhash
    } else { // yes
        MiniRecord mr(prg_id, path, knode_id,
            strand); // create a new MiniRecord from this minimizer kmer
        if (find(it->second->begin(), it->second->end(), mr)
            == it->second->end()) { // checks if this minimizer record is already in the
                                    // vector of this minimizer
            it->second->push_back(mr); // no, add it
        }
    }
}

void Index::clear()
{
    for (auto it = minhash.begin(); it != minhash.end();) {
        delete it->second;
        it = minhash.erase(it);
    }
}

void Index::save(const fs::path& prgfile, uint32_t w, uint32_t k)
{
    const fs::path filepath { prgfile.string() + ".k" + std::to_string(k) + ".w"
        + std::to_string(w) + ".idx" };
    save(filepath);
}

void Index::save(const fs::path& indexfile)
{
    BOOST_LOG_TRIVIAL(debug) << "Saving index to " << indexfile;
    fs::ofstream handle;
    handle.open(indexfile);

    handle << minhash.size() << std::endl;

    for (auto& it : minhash) {
        handle << it.first << "\t" << it.second->size();
        for (uint32_t j = 0; j != it.second->size(); ++j) {
            handle << "\t" << (*(it.second))[j];
        }
        handle << std::endl;
    }
    handle.close();
    BOOST_LOG_TRIVIAL(debug) << "Finished saving " << minhash.size()
                             << " entries to file";
}

void Index::load(fs::path prgfile, uint32_t w, uint32_t k)
{
    const auto ext { ".k" + std::to_string(k) + ".w" + std::to_string(w) + ".idx" };
    prgfile += ext;
    load(prgfile);
}

void Index::load(const fs::path& indexfile)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading index";
    BOOST_LOG_TRIVIAL(debug) << "File is " << indexfile;
    uint64_t key;
    size_t size;
    int c;
    MiniRecord mr;
    bool first = true;

    fs::ifstream myfile(indexfile);
    if (myfile.is_open()) {
        while (myfile.good()) {
            c = myfile.peek();
            if (isdigit(c) and first) {
                myfile >> size;
                minhash.reserve(minhash.size() + size);
                first = false;
                myfile.ignore(1, '\t');
            } else if (isdigit(c) and !first) {
                myfile >> key;
                myfile.ignore(1, '\t');
                myfile >> size;
                auto* vmr = new std::vector<MiniRecord>;
                if (minhash.find(key) != minhash.end()) {
                    vmr = minhash[key];
                    vmr->reserve(vmr->size() + size);
                } else {
                    vmr->reserve(size);
                    minhash[key] = vmr;
                }
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
        fatal_error("Unable to open index file ", indexfile,
            ". Does it exist? Have you run pandora index?");
    }

    if (minhash.size() <= 1) {
        BOOST_LOG_TRIVIAL(debug)
            << "Was this file empty?! Index now contains a trivial " << minhash.size()
            << " entries";
    } else {
        BOOST_LOG_TRIVIAL(debug) << "Finished loading file. Index now contains "
                                 << minhash.size() << " entries";
    }
}

bool Index::operator==(const Index& other) const
{
    if (this->minhash.size() != other.minhash.size()) {
        return false;
    }

    for (const auto& kmer : this->minhash) {
        const auto it = other.minhash.find(kmer.first);
        if (it == other.minhash.end()) {
            return false;
        }
        const auto& qvecp = it->second;
        for (const auto& record : *this->minhash.at(kmer.first)) {
            if (std::find(qvecp->begin(), qvecp->end(), record) == qvecp->end()) {
                return false;
            }
        }
    }

    for (const auto& kmer : other.minhash) {
        const auto it = this->minhash.find(kmer.first);
        if (it == this->minhash.end()) {
            return false;
        }
        const auto& qvecp = it->second;
        for (const auto& record : *other.minhash.at(kmer.first)) {
            if (std::find(qvecp->begin(), qvecp->end(), record) == qvecp->end()) {
                return false;
            }
        }
    }
    return true;
}

bool Index::operator!=(const Index& other) const { return !(*this == other); }

void index_prgs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
    std::shared_ptr<Index>& index, const uint32_t w, const uint32_t k,
    const fs::path& outdir, uint32_t threads)
{
    BOOST_LOG_TRIVIAL(debug) << "Index PRGs";
    if (prgs.empty())
        return;

    // first reserve an estimated index size
    uint32_t r = 0;
    for (uint32_t i = 0; i != prgs.size(); ++i) {
        r += prgs[i]->seq.length();
    }
    index->minhash.reserve(r);

    // create the dirs for the index
    const int nbOfGFAsPerDir = 4000;
    for (uint32_t i = 0; i <= prgs.size() / nbOfGFAsPerDir; ++i)
        fs::create_directories(outdir / int_to_string(i + 1));

    // now fill index
    std::atomic_uint32_t nbOfPRGsDone { 0 };
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint32_t i = 0; i < prgs.size(); ++i) { // for each prg
        uint32_t dir = i / nbOfGFAsPerDir + 1;
        prgs[i]->minimizer_sketch(
            index, w, k, (((double)(nbOfPRGsDone.load())) / prgs.size()) * 100);
        const auto gfa_file { outdir / int_to_string(dir)
            / (prgs[i]->name + ".k" + std::to_string(k) + ".w" + std::to_string(w)
                + ".gfa") };
        prgs[i]->kmer_prg.save(gfa_file);

        ++nbOfPRGsDone;
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished adding " << prgs.size() << " LocalPRGs";
    BOOST_LOG_TRIVIAL(debug) << "Number of keys in Index: " << index->minhash.size();
}
