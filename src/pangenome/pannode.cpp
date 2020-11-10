#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "pangenome/panread.h"
#include "minihit.h"
#include "utils.h"
#include "localPRG.h"
#include "OptionsAggregator.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace pangenome;

// constructors
pangenome::Node::Node(const std::shared_ptr<LocalPRG>& prg)
    : Node(prg, prg->id)
{
}
pangenome::Node::Node(const std::shared_ptr<LocalPRG>& prg, uint32_t node_id,
    uint32_t total_number_samples // total number of samples that we have in this node
    )
    : prg(prg)
    , prg_id(prg->id)
    , node_id(node_id)
    , name(prg->name)
    , covg(0)
    , kmer_prg_with_coverage(
          const_cast<KmerGraph*>(
              &prg->kmer_prg), // TODO: this is very dangerous,
                               // KmerGraphWithCoverage::kmer_prg must be made const
          total_number_samples)
{
}

/*// copy constructor
Node::Node(const Node& other)
{
    prg_id = other.prg_id;
    node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    reads = other.reads;
    samples = other.samples;
}

// Assignment operator
Node& Node::operator=(const Node& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    prg_id = other.prg_id;
    node_id = other.node_id;
    name = other.name;
    covg = other.covg;
    kmer_prg = other.kmer_prg;
    reads = other.reads;
    samples = other.samples;

    return *this;
}*/

void pangenome::Node::remove_read(ReadPtr r)
{
    // removes single copy of read
    auto it = find(reads.begin(), reads.end(), r);
    if (it != reads.end()) {
        covg -= 1;
        reads.erase(it);
        // it = find(reads.begin(), reads.end(), r);
    }
}

std::string pangenome::Node::get_name() const
{
    if (prg_id != node_id) {
        return name + "." + std::to_string(node_id);
    } else {
        return name;
    }
}

void pangenome::Node::add_path(
    const std::vector<KmerNodePtr>& kmp, const uint32_t& sample_id)
{
    for (uint32_t i = 0; i != kmp.size(); ++i) {
        assert(kmp[i]->id < kmer_prg_with_coverage.kmer_prg->nodes.size()
            and kmer_prg_with_coverage.kmer_prg->nodes[kmp[i]->id] != nullptr);
        kmer_prg_with_coverage.increment_forward_covg(kmp[i]->id, sample_id);
        kmer_prg_with_coverage.increment_reverse_covg(kmp[i]->id, sample_id);
    }
}

void pangenome::Node::get_read_overlap_coordinates(
    std::vector<std::vector<uint32_t>>& read_overlap_coordinates)
{
    read_overlap_coordinates.reserve(reads.size());
    std::vector<uint32_t> coordinate;

    auto read_count = 0;
    for (const auto& read_ptr : reads) {
        read_count++;
        auto hits = read_ptr->get_hits_as_unordered_map();
        if (hits.at(prg_id).size() < 2)
            continue;

        auto hit_ptr_iter = hits.at(prg_id).begin();
        uint32_t start = (*hit_ptr_iter)->get_read_start_position();
        uint32_t end = 0;
        for (const auto& hit_ptr : hits.at(prg_id)) {
            start = std::min(start, hit_ptr->get_read_start_position());
            end = std::max(end,
                hit_ptr->get_read_start_position() + hit_ptr->get_prg_path().length());
        }

        assert(end > start
            or assert_msg("Error finding the read overlap coordinates for node "
                << name << " and read " << read_ptr->id << " (the " << read_count
                << "th on this node)" << std::endl
                << "Found end " << end << " after found start " << start));
        coordinate = { read_ptr->id, start, end, (*hit_ptr_iter)->is_forward() };
        read_overlap_coordinates.push_back(coordinate);
    }

    if (not read_overlap_coordinates.empty()) {
        sort(read_overlap_coordinates.begin(), read_overlap_coordinates.end(),
            [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                for (uint32_t i = 0; i < a.size(); ++i) {
                    if (a[i] != b[i]) {
                        return a[i] < b[i];
                    }
                }
                return false;
            });
    }
}

void pangenome::Node::construct_multisample_vcf(VCF& master_vcf,
    const std::vector<LocalNodePtr>& vcf_reference_path,
    const std::shared_ptr<LocalPRG>& prg, const uint32_t w)
{
    // create a vcf with respect to this ref
    VCF vcf(master_vcf.genotyping_options);
    prg->build_vcf_from_reference_path(vcf, vcf_reference_path);
    vcf.add_samples(master_vcf.samples);

    BOOST_LOG_TRIVIAL(debug) << "Initial build:\n"
                             << vcf.to_string(true, false) << std::endl;

    for (const auto& sample : samples) {
        uint32_t count = 0;
        for (const auto& sample_kmer_path : sample->paths[prg_id]) {
            const auto sample_local_path
                = prg->localnode_path_from_kmernode_path(sample_kmer_path, w);
            if (count == 0) {
                prg->add_new_records_and_genotype_to_vcf_using_max_likelihood_path_of_the_sample(
                    vcf, vcf_reference_path, sample_local_path, sample->name);
                BOOST_LOG_TRIVIAL(debug) << "With sample added:\n"
                                         << vcf.to_string(true, false);
                prg->add_sample_covgs_to_vcf(vcf, kmer_prg_with_coverage,
                    vcf_reference_path, sample->name, sample->sample_id);
                BOOST_LOG_TRIVIAL(debug) << "With sample coverages added:\n"
                                         << vcf.to_string(true, false);
            } else {
                auto path_specific_sample_name = sample->name + std::to_string(count);
                prg->add_new_records_and_genotype_to_vcf_using_max_likelihood_path_of_the_sample(
                    vcf, vcf_reference_path, sample_local_path,
                    path_specific_sample_name);
                prg->add_sample_covgs_to_vcf(vcf, kmer_prg_with_coverage,
                    vcf_reference_path, path_specific_sample_name, sample->sample_id);
            }
            count++;
        }
    }
    vcf = vcf.merge_multi_allelic();
    BOOST_LOG_TRIVIAL(debug) << "After merging alleles:\n"
                             << vcf.to_string(true, false);
    vcf = vcf.correct_dot_alleles(
        prg->string_along_path(vcf_reference_path), prg->name);
    BOOST_LOG_TRIVIAL(debug) << "After fixing dot alleles:\n"
                             << vcf.to_string(true, false);
    master_vcf.append_vcf(vcf);
}

bool pangenome::Node::operator==(const Node& y) const { return (node_id == y.node_id); }

bool pangenome::Node::operator!=(const Node& y) const { return (node_id != y.node_id); }

bool pangenome::Node::operator<(const Node& y) const { return (node_id < y.node_id); }

std::ostream& pangenome::operator<<(std::ostream& out, pangenome::Node const& n)
{
    out << n.node_id << "," << n.prg_id << " covg: " << n.covg;
    return out;
}

std::set<ReadCoordinate> pangenome::Node::get_read_overlap_coordinates(
    const prg::Path& local_path, const uint32_t& min_number_hits)
{
    std::set<ReadCoordinate> read_overlap_coordinates;

    for (const auto& current_read : this->reads) {
        auto hits = current_read->get_hits_as_unordered_map();
        const auto read_hits_inside_path { find_hits_inside_path(
            hits.at(this->prg_id), local_path) };

        if (read_hits_inside_path.size() < min_number_hits) {
            continue;
        }

        const auto read_hits_iter { read_hits_inside_path.cbegin() };
        uint32_t start { (*read_hits_iter)->get_read_start_position() };
        uint32_t end { 0 };

        for (const auto& read_hit : read_hits_inside_path) {
            start = std::min(start, read_hit->get_read_start_position());
            end = std::max(end,
                read_hit->get_read_start_position()
                    + read_hit->get_prg_path().length());
        }

        assert(end > start);

        read_overlap_coordinates.emplace(
            current_read->id, start, end, (*read_hits_iter)->is_forward());
    }
    return read_overlap_coordinates;
}