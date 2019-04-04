#ifndef PANDORA_CANDIDATE_REGION_H
#define PANDORA_CANDIDATE_REGION_H

#include <vector>
#include <algorithm>
#include <set>
#include "interval.h"
#include "denovo_utils.h"
#include "fastaq.h"
#include "fastaq_handler.h"
#include "localPRG.h"


using CandidateRegionIdentifier = std::tuple<Interval, std::string>;


template<>
struct std::hash<CandidateRegionIdentifier> {
    size_t operator()(const CandidateRegionIdentifier &id) const;
};


struct TmpPanNode {
    const PanNodePtr &pangraph_node;
    const std::shared_ptr<LocalPRG> &local_prg;
    const std::vector<KmerNodePtr> &kmer_node_max_likelihood_path;
    const std::vector<LocalNodePtr> &local_node_max_likelihood_path;

    TmpPanNode(const PanNodePtr &pangraph_node, const std::shared_ptr<LocalPRG> &local_prg,
               const std::vector<KmerNodePtr> &kmer_node_max_likelihood_path,
               const std::vector<LocalNodePtr> &local_node_max_likelihood_path);
};


using DenovoPaths = std::vector<std::string>;
using ReadPileup = std::vector<std::string>;
using ReadCoordinates = std::set<ReadCoordinate>;


class CandidateRegion {
public:
    ReadCoordinates read_coordinates;
    std::string max_likelihood_sequence;
    std::string left_flanking_sequence;
    std::string right_flanking_sequence;
    ReadPileup pileup;
    fs::path filename;
    DenovoPaths denovo_paths;

    CandidateRegion(const Interval &interval, std::string name);

    CandidateRegion(const Interval &interval, std::string name, const uint_least16_t &interval_padding);

    bool operator==(const CandidateRegion &rhs) const;

    bool operator!=(const CandidateRegion &rhs) const;

    Interval get_interval() const;

    const std::string &get_name() const;

    CandidateRegionIdentifier get_id() const;

    std::string get_max_likelihood_sequence_with_flanks() const;

    void generate_read_pileup(const fs::path &reads_filepath);

    void write_denovo_paths_to_file(const fs::path &output_directory);

private:
    const Interval interval;
    const std::string name;
    const uint_least16_t interval_padding;

    void initialise_filename();

    Fastaq generate_fasta_for_denovo_paths();
};


std::vector<Interval> identify_low_coverage_intervals(const std::vector<uint32_t> &covg_at_each_position,
                                                      const uint32_t &min_required_covg = 2,
                                                      const uint32_t &min_length = 5);

using CandidateRegions = std::unordered_map<CandidateRegionIdentifier, CandidateRegion>;

CandidateRegions find_candidate_regions_for_pan_node(const TmpPanNode &pangraph_node_components,
                                                     const uint_least16_t &candidate_region_interval_padding = 0);

#endif //PANDORA_CANDIDATE_REGION_H
