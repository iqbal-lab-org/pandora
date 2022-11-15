#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <atomic>

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
    for (auto hash_and_minirecords : minhash) {
        delete hash_and_minirecords.second;
    }
    minhash.clear();
}

void Index::save_minhash()
{
    BOOST_LOG_TRIVIAL(debug) << "Saving minhash...";
    index_archive.prepare_new_entry("_minhash");

    {
        std::stringstream handle;
        handle << minhash.size() << std::endl;
        index_archive.write_data(handle.str(), false);
    }


    for (auto& it : minhash) {
        std::stringstream handle;
        handle << it.first << "\t" << it.second->size();
        for (const MiniRecord &mini_record : *it.second) {
            handle << "\t" << mini_record;
        }
        handle << std::endl;
        index_archive.write_data(handle.str(), false);
    }

    BOOST_LOG_TRIVIAL(debug) << "Saving minhash - done!";
}

Index Index::load(const fs::path& indexfile)
{
    /*
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
     */
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

void Index::index_prgs(LocalPRGReaderGeneratorIterator &prg_it, uintmax_t estimated_index_size,
    const uint32_t indexing_upper_bound, const uint32_t threads)
{
    BOOST_LOG_TRIVIAL(debug) << "Index PRGs";
    index_archive.prepare_new_entry("kmer_prgs.mgfas");

    this->minhash.clear();
    this->minhash.reserve(estimated_index_size);

    // now fill index
    std::atomic_uint32_t nb_of_prgs_done { 0 };
#pragma omp parallel num_threads(threads)
    while (true) {
        std::shared_ptr<LocalPRG> local_prg(nullptr);
#pragma omp critical(prg_it)
        {
            local_prg = *prg_it;
            ++prg_it;
        }

        const bool no_more_local_prgs_to_process = local_prg == nullptr;
        if (no_more_local_prgs_to_process) {
            break;
        }

        try {
            local_prg->minimizer_sketch(this, w, k, indexing_upper_bound, -1);
            std::string prg_as_gfa = local_prg->kmer_prg.to_gfa(local_prg);
            index_archive.write_data(prg_as_gfa, true);
        } catch (const IndexingLimitReached &error) {
            BOOST_LOG_TRIVIAL(warning) << error.what();
        }

        ++nb_of_prgs_done;
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished adding " << nb_of_prgs_done << " LocalPRGs";
    BOOST_LOG_TRIVIAL(debug) << "Number of keys in Index: " << this->minhash.size();

    this->save_minhash();
}

void Index::index_prgs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
    const uint32_t indexing_upper_bound, const uint32_t threads)
{
    BOOST_LOG_TRIVIAL(debug) << "Index PRGs";
    index_archive.prepare_new_entry("kmer_prgs.mgfas");

    this->minhash.clear();

    if (prgs.empty())
        return;

    // first reserve an estimated index size
    uint32_t estimated_index_size = 0;
    for (const auto &prg : prgs) {
        estimated_index_size += prg->seq.length();
    }
    this->minhash.reserve(estimated_index_size);


    // now fill index
    std::atomic_uint32_t nb_of_prgs_done { 0 };
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1) default(shared)
    for (uint32_t i = 0; i < prgs.size(); ++i) { // for each prg
        std::shared_ptr<LocalPRG>& local_prg(prgs[i]);
        try {
            local_prg->minimizer_sketch(this, w, k, indexing_upper_bound,
                (((double)(nb_of_prgs_done.load())) / prgs.size()) * 100);
            std::string prg_as_gfa = local_prg->kmer_prg.to_gfa(local_prg);
            index_archive.write_data(prg_as_gfa, true);
        } catch (const IndexingLimitReached &error) {
            BOOST_LOG_TRIVIAL(warning) << error.what();
        }

        // saves some RAM
        local_prg.reset();

        ++nb_of_prgs_done;
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished adding " << prgs.size() << " LocalPRGs";
    BOOST_LOG_TRIVIAL(debug) << "Number of keys in Index: " << this->minhash.size();

    this->save_minhash();
}
