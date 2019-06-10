#ifndef __KMERNODE_H_INCLUDED__   // if kmernode.h hasn't been included yet...
#define __KMERNODE_H_INCLUDED__

class KmerNode;

#include <cstring>
#include <cstdint>
#include <vector>
#include <ostream>
#include <memory>
#include <algorithm>

#include "prg/path.h"
#include "pangenome/ns.cpp"


typedef std::shared_ptr<KmerNode> KmerNodePtr;
typedef std::weak_ptr<KmerNode> WeakKmerNodePtr;

class KmerNode { //represent a kmer-minimizer in the KmerGraph
private:
    std::vector<WeakKmerNodePtr>::const_iterator findNodePtrInNodesVector(const std::vector<WeakKmerNodePtr> &nodesVector, const KmerNodePtr &rhs) const {
        find_if(nodesVector.begin(), nodesVector.end(), [&rhs](const WeakKmerNodePtr &lhs) {
            return lhs.lock() == rhs;
        });
    }

    std::vector<WeakKmerNodePtr>::const_iterator findNodeInNodesVector(const std::vector<WeakKmerNodePtr> &nodesVector, const KmerNode &rhs) const {
        find_if(nodesVector.begin(), nodesVector.end(), [&rhs](const WeakKmerNodePtr &lhs) {
            return *(lhs.lock()) == rhs;
        });
    }



public:
    //attributes (TODO: protect these? only this class should operate in these attributes, move logic that change them to here?)
    uint32_t id;
    prg::Path path; //the path of the kmer in the LocalPRG
    std::vector<WeakKmerNodePtr> outNodes; // representing edges from this node to the nodes in the vector
    std::vector<WeakKmerNodePtr> inNodes; // representing edges from other nodes to this node
    uint64_t khash; //the kmer hash value
    uint8_t num_AT; // the number of As and Ts in this kmer

    std::vector<WeakKmerNodePtr>::const_iterator findNodePtrInOutNodes(const KmerNodePtr &rhs) const {
        return findNodePtrInNodesVector(outNodes, rhs);
    }
    std::vector<WeakKmerNodePtr>::const_iterator findNodePtrInInNodes(const KmerNodePtr &rhs) const {
        return findNodePtrInNodesVector(inNodes, rhs);
    }
    std::vector<WeakKmerNodePtr>::const_iterator findNodeInOutNodes(const KmerNode &rhs) const {
        return findNodeInNodesVector(outNodes, rhs);
    }

    std::vector<WeakKmerNodePtr>::const_iterator findNodeInInNodes(const KmerNode &rhs) const {
        return findNodeInNodesVector(inNodes, rhs);
    }

    //constructors and assignment operators
    KmerNode(uint32_t, const prg::Path &);
    KmerNode(const KmerNode &);
    KmerNode &operator=(const KmerNode &);
    bool operator==(const KmerNode &y) const;


    //friends
    friend std::ostream &operator<<(std::ostream &out, const KmerNode &n);
    friend class KmerGraph;
    friend struct condition;
    friend struct pCompKmerNode;
    friend class LocalPRG;
    friend class pangenome::Graph;
    friend class pangenome::Node;
    friend int pandora_check_kmergraph(int argc, char *argv[]);

    /*
    friend void estimate_parameters(pangenome::Graph *, const std::string &, const uint32_t, float &, const uint32_t);
     */

};

#endif
