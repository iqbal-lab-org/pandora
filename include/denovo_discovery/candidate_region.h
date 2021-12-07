#ifndef PANDORA_CANDIDATE_REGION_H
#define PANDORA_CANDIDATE_REGION_H

#include "denovo_utils.h"
#include "fastaq.h"
#include "fastaq_handler.h"
#include "interval.h"
#include "localPRG.h"
#include <algorithm>
#include <set>
#include <vector>
#include <sstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#ifndef NO_OPENMP
#include <omp.h>
#endif

using CandidateRegionIdentifier = std::tuple<Interval, std::string>;

template <> struct std::hash<CandidateRegionIdentifier> {
    size_t operator()(const CandidateRegionIdentifier& id) const;
};

struct TmpPanNode {
    const PanNodePtr& pangraph_node;
    const std::shared_ptr<LocalPRG>& local_prg;
    const std::vector<KmerNodePtr>& kmer_node_max_likelihood_path;
    const std::vector<LocalNodePtr>& local_node_max_likelihood_path;
    const std::string& local_node_max_likelihood_seq;

    TmpPanNode(const PanNodePtr& pangraph_node,
        const std::shared_ptr<LocalPRG>& local_prg,
        const std::vector<KmerNodePtr>& kmer_node_max_likelihood_path,
        const std::vector<LocalNodePtr>& local_node_max_likelihood_path,
        const std::string& local_node_max_likelihood_seq);
};

using DenovoPaths = std::vector<std::string>;
using ReadPileup = std::vector<std::string>;
using ReadCoordinates = std::set<ReadCoordinate>;

class CandidateRegionWriteBuffer;

class CandidateRegion {
public:
    ReadCoordinates read_coordinates;
    std::string max_likelihood_sequence;

    const std::vector<LocalNodePtr> local_node_max_likelihood_path;
    // NB local_node_max_likelihood_sequence differs from max_likelihood_sequence in the
    // sense that max_likelihood_sequence is the ML sequence of the candidate region
    // slice only
    const std::string local_node_max_likelihood_sequence;

    std::string left_flanking_sequence;
    std::string right_flanking_sequence;
    ReadPileup pileup;
    fs::path filename;
    DenovoPaths denovo_paths;

    CandidateRegion(const Interval& interval, std::string name);

    CandidateRegion(const Interval& interval, std::string name,
        const uint_least16_t& interval_padding);

    CandidateRegion(const Interval& interval, std::string name,
        const uint_least16_t& interval_padding,
        const std::vector<LocalNodePtr>& local_node_max_likelihood_path,
        const std::string& local_node_max_likelihood_sequence);

    virtual ~CandidateRegion();

    bool operator==(const CandidateRegion& rhs) const;

    bool operator!=(const CandidateRegion& rhs) const;

    Interval get_interval() const;

    const std::string& get_name() const;

    CandidateRegionIdentifier get_id() const;

    std::string get_max_likelihood_sequence_with_flanks() const;

    void add_pileup_entry(
        const std::string& read, const ReadCoordinate& read_coordinate);

    virtual std::vector<std::string> get_variants(const string& denovo_sequence) const;

    void write_denovo_paths_to_buffer(CandidateRegionWriteBuffer& buffer);

protected:
    const Interval interval;
    const std::string name;
    const uint_least16_t interval_padding;

#ifndef NO_OPENMP
    // TODO: check if we might have issues here with copy constructor - I think not:
    // OpenMP locks work on the address of the lock I think
    omp_lock_t add_pileup_entry_lock; // synchronizes multithreaded access to the
                                      // add_pileup_entry() method
#endif

    void init();
    void initialise_filename();
};

using CandidateRegions = std::unordered_map<CandidateRegionIdentifier, CandidateRegion>;
using ReadId = uint32_t;
using PileupConstructionMap
    = std::map<ReadId, std::vector<std::pair<CandidateRegion*, const ReadCoordinate*>>>;

class Discover {
private:
    const uint32_t min_required_covg;
    const uint32_t min_candidate_len;
    const uint32_t max_candidate_len;
    const uint16_t candidate_padding;
    const uint32_t merge_dist;

public:
    Discover(uint32_t min_required_covg, uint32_t min_candidate_len,
        uint32_t max_candidate_len, uint16_t candidate_padding, uint32_t merge_dist);

    std::vector<Interval> identify_low_coverage_intervals(
        const std::vector<uint32_t>& covg_at_each_position);

    CandidateRegions find_candidate_regions_for_pan_node(
        const TmpPanNode& pangraph_node_components);

    PileupConstructionMap pileup_construction_map(CandidateRegions& candidate_regions);

    void load_candidate_region_pileups(const fs::path& reads_filepath,
        const CandidateRegions& candidate_regions,
        const PileupConstructionMap& pileup_construction_map, uint32_t threads = 1);
};

class CandidateRegionWriteBuffer {
protected:
    std::string sample_name;
    std::map<std::string, std::string> locus_name_to_ML_path;
    std::map<std::string, std::vector<std::string>> locus_name_to_variants;

    template <typename OFSTREAM_TYPE>
    void write_to_file_core(OFSTREAM_TYPE& output_filehandler) const
    {
        output_filehandler << "Sample " << sample_name << std::endl;
        output_filehandler << locus_name_to_ML_path.size()
                           << " loci with denovo variants" << std::endl;
        for (const auto& locus_name_and_ML_path : locus_name_to_ML_path) {
            const auto& locus_name = locus_name_and_ML_path.first;
            output_filehandler << locus_name << std::endl;
            output_filehandler << locus_name_and_ML_path.second << std::endl;
            output_filehandler << locus_name_to_variants.at(locus_name).size()
                               << " denovo variants for this locus" << std::endl;
            for (const auto& variant : locus_name_to_variants.at(locus_name)) {
                output_filehandler << variant << std::endl;
            }
        }
    }

    friend class boost::serialization::access;
    // Note: trivial method, not tested
    template <class Archive> void serialize(Archive& ar, const unsigned int version)
    {
        ar& sample_name;
        ar& locus_name_to_ML_path;
        ar& locus_name_to_variants;
    }

public:
    CandidateRegionWriteBuffer() { } // required by Boost serialization

    CandidateRegionWriteBuffer(const std::string& sample_name)
        : sample_name(sample_name)
    {
    }

    virtual ~CandidateRegionWriteBuffer() = default;
    CandidateRegionWriteBuffer(const CandidateRegionWriteBuffer& other) = default;
    CandidateRegionWriteBuffer(CandidateRegionWriteBuffer&& other) = default;
    CandidateRegionWriteBuffer& operator=(const CandidateRegionWriteBuffer& other)
        = default;
    CandidateRegionWriteBuffer& operator=(CandidateRegionWriteBuffer&& other) = default;
    bool operator==(const CandidateRegionWriteBuffer& other) const
    {
        return this->sample_name == other.sample_name
            && this->locus_name_to_ML_path == other.locus_name_to_ML_path
            && this->locus_name_to_variants == other.locus_name_to_variants;
    }
    bool operator!=(const CandidateRegionWriteBuffer& other) const
    {
        return !(*this == other);
    }

    virtual void add_new_variant(const std::string& locus_name,
        const std::string& ML_path, const std::string& variant);

    void write_to_file(const fs::path& output_file) const;

    const std::string& get_sample_name() const { return sample_name; }

    const std::map<std::string, std::string>& get_locus_name_to_ML_path() const
    {
        return locus_name_to_ML_path;
    }

    const std::map<std::string, std::vector<std::string>>
    get_locus_name_to_variants() const
    {
        return locus_name_to_variants;
    }

    void merge(const CandidateRegionWriteBuffer& other);
};

#endif // PANDORA_CANDIDATE_REGION_H
