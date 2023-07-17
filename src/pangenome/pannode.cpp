#include <iostream>
#include <algorithm>
#include <boost/log/trivial.hpp>
#include "pangenome/pannode.h"
#include "pangenome/pansample.h"
#include "minihit.h"
#include "localPRG.h"
#include "OptionsAggregator.h"
#include "denovo_discovery/denovo_utils.h"

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
        const bool kmer_node_is_valid
            = (kmp[i]->id < kmer_prg_with_coverage.kmer_prg->nodes.size())
            and (kmer_prg_with_coverage.kmer_prg->nodes[kmp[i]->id] != nullptr);
        if (!kmer_node_is_valid) {
            fatal_error(
                "When adding a path to a Pangraph Node, a kmer node is not valid");
        }
        kmer_prg_with_coverage.increment_forward_covg(kmp[i]->id, sample_id);
        kmer_prg_with_coverage.increment_reverse_covg(kmp[i]->id, sample_id);
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
