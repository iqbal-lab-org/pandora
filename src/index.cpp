#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <atomic>

#include <boost/log/trivial.hpp>

#include "minirecord.h"
#include "index.h"
#include "localPRG.h"

std::vector<std::string> get_prg_names(
    const fs::path &prg_filepath, uintmax_t estimated_index_size) {
    std::vector<std::string> prg_names;

    // we estimate the number of PRGs thinking each PRG has the size of ~1000 chars/bytes
    // not the best but not the worst estimation
    prg_names.reserve(estimated_index_size/1000);

    std::ifstream prg_ifstream;
    open_file_for_reading(prg_filepath.string(), prg_ifstream);

    std::string line;
    while (std::getline(prg_ifstream, line))
    {
        if (!line.empty()) {
            const bool is_header = line[0] == '>';
            if (is_header) {
                const std::string prg_name = line.substr(1);  // remove '>'
                prg_names.push_back(prg_name);
            }
        }
    }

    prg_ifstream.close();

    return prg_names;
}

void Index::build_index_on_disk(const uint32_t w, const uint32_t k, const fs::path &prg_filepath,
    const fs::path &out_filepath, const uint32_t indexing_upper_bound,
    const uint32_t threads) {
    ZipFileWriter index_archive(out_filepath);

    uintmax_t estimated_index_size = get_number_of_bytes_in_file(prg_filepath);
    std::vector<std::string> prg_names = get_prg_names(prg_filepath, estimated_index_size);

    Index index(w, k, std::move(prg_names));

    // load PRGs from file lazily (using generators)
    LocalPRGGeneratorAndIterator prg_generator_and_it = LocalPRGReader::read_prg_file_as_generator(prg_filepath);

    index.index_prgs(index_archive, prg_generator_and_it.second,
        estimated_index_size, indexing_upper_bound, threads);
}

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

void Index::save_minhash(ZipFileWriter &index_archive) const
{
    BOOST_LOG_TRIVIAL(debug) << "Saving minhash...";
    index_archive.prepare_new_entry("_minhash");

    {
        std::stringstream handle;
        handle << minhash.size() << std::endl;
        index_archive.write_data(handle.str());
    }


    for (const auto& it : minhash) {
        std::stringstream handle;
        handle << it.first << "\t" << it.second->size();
        for (const MiniRecord &mini_record : *it.second) {
            handle << "\t" << mini_record;
        }
        handle << std::endl;
        index_archive.write_data(handle.str());
    }

    BOOST_LOG_TRIVIAL(debug) << "Saving minhash - done!";
}

void Index::load_minhash(ZipFileReader &zip_file) {
    uint64_t key;
    size_t size;
    int c;
    MiniRecord mr;
    bool first = true;

    ZipInstreamBuffer zip_buffer(zip_file.archive, "_minhash");
    std::istream zip_in(&zip_buffer);

    while (zip_in.good()) {
        c = zip_in.peek();
        if (isdigit(c) and first) {
            zip_in >> size;
            minhash.reserve(minhash.size() + size);
            first = false;
            zip_in.ignore(1, '\t');
        } else if (isdigit(c) and !first) {
            zip_in >> key;
            zip_in.ignore(1, '\t');
            zip_in >> size;
            auto* vmr = new std::vector<MiniRecord>;
            if (minhash.find(key) != minhash.end()) {
                vmr = minhash[key];
                vmr->reserve(vmr->size() + size);
            } else {
                vmr->reserve(size);
                minhash[key] = vmr;
            }
            zip_in.ignore(1, '\t');
        } else if (c == EOF) {
            break;
        } else {
            zip_in >> mr;
            minhash[key]->push_back(mr);
            zip_in.ignore(1, '\t');
        }
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

Index Index::load(const fs::path& indexfile, std::vector<std::shared_ptr<LocalPRG>>& prgs)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading index";
    BOOST_LOG_TRIVIAL(debug) << "File is " << indexfile;

    ZipFileReader zip_file(indexfile);

    BOOST_LOG_TRIVIAL(debug) << "Loading metadata";
    auto w_and_k = zip_file.read_w_and_k();

    BOOST_LOG_TRIVIAL(debug) << "Loading prg names";
    auto prg_names = zip_file.read_prg_names();

    BOOST_LOG_TRIVIAL(debug) << "Loading prg min path lengths";
    auto prg_min_path_lengths = zip_file.read_prg_min_path_lengths();
    Index index(w_and_k.first, w_and_k.second, std::move(prg_names),
        std::move(prg_min_path_lengths));

    BOOST_LOG_TRIVIAL(debug) << "Loading minhash";
    index.load_minhash(zip_file);

    BOOST_LOG_TRIVIAL(debug) << "Loading kmer prgs";
    for (const auto& prg : prgs) {
        const auto filename { prg->name + ".gfa" };
        std::string gfa_as_str = zip_file.read_full_text_file_as_single_string(filename);
        std::stringstream gfa_ss;
        gfa_ss << gfa_as_str;
        prg->kmer_prg.load(gfa_ss);
    }

    return index;
}

bool Index::operator==(const Index& other) const
{
    if (this->w != other.w) {
        return false;
    }

    if (this->k != other.k) {
        return false;
    }

    // TODO: check prg_names and prg_min_path_lengths?

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

std::unordered_map<std::string, uint32_t> Index::get_prg_names_to_prg_index() const {
    std::unordered_map<std::string, uint32_t> prg_names_to_prg_index;
    prg_names_to_prg_index.reserve(prg_names.size() * 2);
    for (uint32_t prg_index=0; prg_index < prg_names.size(); prg_index++) {
        prg_names_to_prg_index[prg_names[prg_index]] = prg_index;
    }
    return prg_names_to_prg_index;
}

void Index::index_prgs(ZipFileWriter &index_archive,
    LocalPRGReaderGeneratorIterator &prg_it, uintmax_t estimated_index_size,
    const uint32_t indexing_upper_bound, const uint32_t threads)
{
    BOOST_LOG_TRIVIAL(debug) << "Index PRGs";

    // minhash init
    this->minhash.reserve(estimated_index_size);

    // builds prg_names_to_prg_index to be used in the main loop
    const std::unordered_map<std::string, uint32_t> prg_names_to_prg_index = get_prg_names_to_prg_index();

    // fill prg_min_path_lengths with 0-ed values that will be set in the main loop
    prg_min_path_lengths.resize(prg_names.size());

    // now fill index
#pragma omp parallel num_threads(threads)
    while (true) {
        std::shared_ptr<LocalPRG> local_prg(nullptr);

        // gets a new PRG
#pragma omp critical(Index__index_prgs__prg_it)
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

            // zip file variables
            std::string prg_as_gfa = local_prg->kmer_prg.to_gfa();
            std::string zip_path = local_prg->name + ".gfa";

            // prg_min_path_lengths variables
            uint32_t prg_index = prg_names_to_prg_index.at(local_prg->name);
            uint32_t min_path_length = local_prg->kmer_prg.min_path_length();

#pragma omp critical(Index__index_prgs__write_gfa_to_zip)
            {
                index_archive.prepare_new_entry(zip_path);
                index_archive.write_data(prg_as_gfa);

            }

#pragma omp critical(Index__index_prgs__prg_min_path_lengths)
            {
                prg_min_path_lengths[prg_index] = min_path_length;
            }
        } catch (const IndexingLimitReached &error) {
            BOOST_LOG_TRIVIAL(warning) << error.what();
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Number of keys in index: " << this->minhash.size();

    // save remaining data
    this->save_minhash(index_archive);
    this->save_prg_names(index_archive);
    this->save_prg_min_path_lengths(index_archive);
    this->save_metadata(index_archive);
}
