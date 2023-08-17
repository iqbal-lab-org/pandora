#ifndef __PANNODE_H_INCLUDED__ // if pannode.h hasn't been included yet...
#define __PANNODE_H_INCLUDED__

#include <string>
#include <cstdint>
#include <unordered_set>
#include <vector>
#include "kmergraph.h"
#include "kmergraphwithcoverage.h"
#include "localPRG.h"
#include "vcf.h"
#include "pansample.h"
#include "OptionsAggregator.h"

class LocalPRG;

class pangenome::Node {
public:
    std::set<SamplePtr, SamplePtrSorterBySampleId> samples;
    const uint32_t prg_id; // corresponding the the LocalPRG id - TODO: this is not
                           // needed - we point to the LocalPRG, which has this info
    const uint32_t
        node_id; // unique node id, so can have multiple copies of a localPRG in graph
    const std::string name; // TODO: this is not needed - we point to the LocalPRG,
                            // which has this info
    mutable uint32_t covg;
    std::shared_ptr<LocalPRG> prg; // TODO: this should be made const
    KmerGraphWithCoverage kmer_prg_with_coverage;

    // main constructor
    Node(const std::shared_ptr<LocalPRG>& prg, uint32_t node_id,
        uint32_t total_number_samples
        = 1 // total number of samples that we have in this node
    );

    // convenience constructors
    Node(const std::shared_ptr<LocalPRG>& prg);

    // Node(const Node&);
    // Node& operator=(const Node&);

    std::string get_name() const;

    void add_path(const std::vector<KmerNodePtr>&, const uint32_t& sample_id);

    void construct_multisample_vcf(VCF& master_vcf,
        const std::vector<LocalNodePtr>& vcf_reference_path,
        const std::shared_ptr<LocalPRG>& prg, const uint32_t w);

    bool operator==(const Node& y) const;

    bool operator!=(const Node& y) const;

    bool operator<(const Node& y) const;

    friend std::ostream& operator<<(std::ostream& out, const Node& n);

    friend class pangenome::Graph;
};

struct EqualComparatorWeakNodePtr {
    bool operator()(
        const pangenome::WeakNodePtr& lhs, const pangenome::WeakNodePtr& rhs)
    {
        return *(lhs.lock()) == *(rhs.lock());
    }
};

#endif
