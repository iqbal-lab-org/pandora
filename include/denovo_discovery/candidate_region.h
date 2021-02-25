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
    const std::string &local_node_max_likelihood_seq;

    TmpPanNode(const PanNodePtr& pangraph_node,
        const std::shared_ptr<LocalPRG>& local_prg,
        const std::vector<KmerNodePtr>& kmer_node_max_likelihood_path,
        const std::vector<LocalNodePtr>& local_node_max_likelihood_path,
        const std::string &local_node_max_likelihood_seq);
};

using DenovoPaths = std::vector<std::string>;
using ReadPileup = std::vector<std::string>;
using ReadCoordinates = std::set<ReadCoordinate>;

class CandidateRegionWriteBuffer;

class CandidateRegion {
public:
    ReadCoordinates read_coordinates;
    std::string max_likelihood_sequence;
    std::vector<LocalNodePtr> local_node_max_likelihood_path;
    std::string local_node_max_likelihood_sequence;
    std::string left_flanking_sequence;
    std::string right_flanking_sequence;
    ReadPileup pileup;
    fs::path filename;
    DenovoPaths denovo_paths;

    CandidateRegion(const Interval& interval, std::string name);

    CandidateRegion(const Interval& interval, std::string name,
        const uint_least16_t& interval_padding,
        const std::vector<LocalNodePtr> &local_node_max_likelihood_path,
        const std::string &local_node_max_likelihood_sequence);

    virtual ~CandidateRegion();

    bool operator==(const CandidateRegion& rhs) const;

    bool operator!=(const CandidateRegion& rhs) const;

    Interval get_interval() const;

    const std::string& get_name() const;

    CandidateRegionIdentifier get_id() const;

    std::string get_max_likelihood_sequence_with_flanks() const;

    void add_pileup_entry(
        const std::string& read, const ReadCoordinate& read_coordinate);

    void write_denovo_paths_to_buffer(CandidateRegionWriteBuffer &buffer);

private:
    const Interval interval;
    const std::string name;
    const uint_least16_t interval_padding;

#ifndef NO_OPENMP
    // TODO: check if we might have issues here with copy constructor - I think not:
    // OpenMP locks work on the address of the lock I think
    omp_lock_t add_pileup_entry_lock; // synchronizes multithreaded access to the
                                      // add_pileup_entry() method
#endif

    void initialise_filename();

    Fastaq generate_fasta_for_denovo_paths();

    template <typename ItType>
    size_t find_length_of_longest_common_prefix(ItType s1_begin, ItType s1_end,
                                                ItType s2_begin, ItType s2_end) const {
        size_t length_of_LCP;
        ItType s1_it, s2_it;
        for (length_of_LCP = 0, s1_it = s1_begin, s2_it = s2_begin;
             s1_it < s1_end && s2_it < s2_end;
             ++s1_it, ++s2_it, ++length_of_LCP) {
            bool bases_are_equal = (*s1_it) == (*s2_it);
            if (!bases_are_equal) {
                break;
            }
        }
        return length_of_LCP;
    }
    std::string get_variants(const string &denovo_path) const;
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
private:
    std::string sample_name;
    std::map<std::string, std::string> locus_name_to_ML_path;
    std::map<std::string, std::vector<std::string>> locus_name_to_variants;

public:
    CandidateRegionWriteBuffer(const std::string &sample_name) :
        sample_name(sample_name){}

    void add_new_variant(const std::string &locus_name,
                         const std::string &ML_path,
                         const std::string &variant);
    void write_to_file(const fs::path& output_file);
};

#endif // PANDORA_CANDIDATE_REGION_H
