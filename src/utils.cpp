#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <cassert>
#include <memory>
//#include "FastxParser.hpp"
#include "utils.h"
#include "seq.h"
#include "pangraph.h"
#include "panread.h"
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

string now() {
    time_t now;
    string dt;

    now = time(nullptr);
    dt = ctime(&now);
    return dt.substr(0, dt.length() - 1) + " ";
}

vector<string> split(const string &query, const string &d) {
    vector<string> v;
    string::size_type k = 0;
    string::size_type j = query.find(d, k);
    while (j != string::npos) {
        if (j > k) {
            v.push_back(query.substr(k, j - k));
        }
        k = j + d.size();
        j = query.find(d, k);
    }
    if (k < query.length()) {
        v.push_back(query.substr(k));
    }
    return v;
}

char complement(char n) {
    switch (n) {
        case 'A':
        case 'a':
            return 'T';
        case 'T':
        case 't':
            return 'A';
        case 'G':
        case 'g':
            return 'C';
        case 'C':
        case 'c':
            return 'G';
    }
    //assert(false);
    return 'N';
}

string rev_complement(string s) {
    transform(begin(s), end(s), begin(s), complement);
    reverse(s.begin(), s.end());
    return s;
}

float lognchoosek2(uint32_t n, uint32_t k1, uint32_t k2) {
    assert(n >= k1 + k2 || assert_msg(
            "Currently the model assumes that the most a given kmer (defined by position) can occur is once per read, i.e. an error somewhere else in the read cannot result in this kmer. If you are getting this message, then you have evidence of violation of this assumption. Either try using a bigger k, or come up with a better model"));
    float total = 0;

    for (uint m = n; m != n - k1 - k2; --m) {
        total += log(m);
    }

    for (uint m = 1; m < k1; ++m) {
        total -= log(m + 1);
    }

    for (uint m = 1; m < k2; ++m) {
        total -= log(m + 1);
    }

    return total;
}

void read_prg_file(vector<LocalPRG *> &prgs, const string &filepath) {
    cout << now() << "Loading PRGs from file " << filepath << endl;

    uint32_t id = 0;
    string name, read, line;
    LocalPRG *s;

    ifstream myfile(filepath);
    if (myfile.is_open()) {
        uint i = 0;
        while (getline(myfile, line).good()) {
            if (line.empty() || line[0] == '>') {
                if (!name.empty() && !read.empty()) {
                    //cout << now() << "Found PRG " << name << endl;
                    s = new LocalPRG(id, name, read);

                    if (s != nullptr) {
                        prgs.push_back(s);
                        id++;
                    } else {
                        cerr << "Failed to make LocalPRG for " << name << endl;
                        exit(1);
                    }
                }
                name.clear();
                read.clear();
                if (!line.empty()) {
                    name = line.substr(1);
                }
            } else {
                read += line;
            }
            i++;
        }
        // and last entry
        if (!name.empty() && !read.empty()) {
            //cout << now() << "Found PRG " << name << endl;
            s = new LocalPRG(id, name, read);
            if (s != nullptr) {
                prgs.push_back(s);
            } else {
                cerr << "Failed to make LocalPRG for " << name << endl;
                exit(1);
            }
        }
        cout << now() << "Number of LocalPRGs read: " << prgs.size() << endl;
        myfile.close();
    } else {
        cerr << "Unable to open PRG file " << filepath << endl;
        exit(EXIT_FAILURE);
    }
}

void load_PRG_kmergraphs(vector<LocalPRG *> &prgs, const uint &w, const uint &k, const string &prgfile) {
    string prefix = "";
    size_t pos = prgfile.find_last_of("/");
    if (pos != std::string::npos) {
        prefix += prgfile.substr(0, pos);
        prefix += "/";
    }
    //cout << "prefix for kmerprgs dir is " << prefix << endl; 
    for (uint i = 0; i != prgs.size(); ++i) {
        prgs[i]->kmer_prg.load(
                prefix + "kmer_prgs/" + prgs[i]->name + ".k" + to_string(k) + ".w" + to_string(w) + ".gfa");
    }
}

//void add_read_hits(const uint32_t id, const string& name, const string& seq, MinimizerHits* hits, Index* idx, const uint32_t w, const uint32_t k)
void add_read_hits(Seq *s, MinimizerHits *hits, Index *idx) {
    //cout << now() << "Search for hits for read " << s->name << " which has sketch size " << s->sketch.size() << " against index of size " << idx->minhash.size() << endl;
    uint32_t hit_count = 0;
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and adds hits to a global MinimizerHits object
    //Seq s(id, name, seq, w, k);
    for (auto it = s->sketch.begin(); it != s->sketch.end(); ++it) {
        if (idx->minhash.find((*it)->kmer) != idx->minhash.end()) {
            for (uint j = 0; j != idx->minhash[(*it)->kmer]->size(); ++j) {
                hits->add_hit(s->id, *it, &(idx->minhash[(*it)->kmer]->operator[](j)));
                hit_count += 1;
            }
            //} else {
            //    cout << "did not find minimizer " << (*it)->kmer << " in index" << endl;
        }
    }
    //hits->sort();
    //cout << now() << "Found " << hit_count << " hits found for read " << s->name << " so size of MinimizerHits is now "
    //     << hits->hits.size() + hits->uhits.size() << endl;
}

void define_clusters(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits, const vector<LocalPRG *> &prgs,
                     const MinimizerHits *minimizer_hits, const int max_diff, const float& scale_cluster_size,
                     const uint min_cluster_size, const uint short_read_length) {
    cout << now() << "Define clusters of hits from the " << minimizer_hits->hits.size() << " hits" << endl;

    if (minimizer_hits->hits.empty()) { return; }

    // A cluster of hits should match same localPRG, each hit not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    auto mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHitPtr, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    uint length_based_threshold;
    for (auto mh_current = ++minimizer_hits->hits.begin();
         mh_current != minimizer_hits->hits.end(); ++mh_current) {
        if ((*mh_current)->read_id != (*mh_previous)->read_id or
            (*mh_current)->prg_id != (*mh_previous)->prg_id or
            (*mh_current)->strand != (*mh_previous)->strand or
            (abs((int) (*mh_current)->read_interval.start - (int) (*mh_previous)->read_interval.start)) > max_diff) {
            // keep clusters which cover at least 10% of the shortest kmer path
            length_based_threshold = min(prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length(),
                                        short_read_length)*scale_cluster_size;
            /*cout << "gene length " << prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()
                 << " short read length: " << short_read_length
                 << " scale cluster size: " << scale_cluster_size
                 << " length based threshold: " << length_based_threshold
                 << " min cluster size: " << min_cluster_size << endl;*/
            if (current_cluster.size() >
                max(length_based_threshold, min_cluster_size)) {
                clusters_of_hits.insert(current_cluster);
	        //cout << "Found cluster of size " << current_cluster.size() << endl;
                /*} else {
                    cout << "rejected hits" << endl;
                    for (set<MinimizerHit*, pComp>::iterator p=current_cluster.begin(); p!=current_cluster.end(); ++p)
                    {
                        cout << **p << endl;
                    }*/
            }
            current_cluster.clear();
        }
        current_cluster.insert(*mh_current);
        mh_previous = mh_current;
    }
    length_based_threshold = min(prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length(),
                                short_read_length)*scale_cluster_size;
    /*cout << "gene length " << prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()
         << " short read length: " << short_read_length
         << " scale cluster size: " << scale_cluster_size
         << " length based threshold: " << length_based_threshold
         << " min cluster size: " << min_cluster_size << endl;*/
    if (current_cluster.size() >
        max(length_based_threshold, min_cluster_size))
    {
        clusters_of_hits.insert(current_cluster);
	    //cout << "Found cluster of size " << current_cluster.size() << endl;
        /*} else {
            cout << "rejected hits" << endl;
            for (set<MinimizerHit*, pComp>::iterator p=current_cluster.begin(); p!=current_cluster.end(); ++p)
            {
                cout << **p << endl;
            }*/
    }

    cout << now() << "Found " << clusters_of_hits.size() << " clusters of hits" << endl;
}

void filter_clusters(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits) {
    // Next order clusters, choose between those that overlap by too much
    cout << now() << "Filter the " << clusters_of_hits.size() << " clusters of hits " << endl;
    if (clusters_of_hits.empty()) { return; }
    // to do this consider pairs of clusters in turn
    auto c_previous = clusters_of_hits.begin();
    /*cout << "first cluster" << endl;
    for (set<MinimizerHit*, pComp>::iterator p=c_previous->begin(); p!=c_previous->end(); ++p)
    {
        cout << **p << endl;
    }*/
    for (auto c_current = ++clusters_of_hits.begin();
         c_current != clusters_of_hits.end(); ++c_current) {
        /*cout << "current cluster" << endl;
        for (set<MinimizerHit*, pComp>::iterator p=c_current->begin(); p!=c_current->end(); ++p)
        {
            cout << **p << endl;
        }*/
        if (((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) && // if on same read and either
            ((((*(*c_current).begin())->prg_id == (*(*c_previous).begin())->prg_id) && // same prg, different strand
              ((*(*c_current).begin())->strand != (*(*c_previous).begin())->strand)) or // or cluster is contained
             ((*--(*c_current).end())->read_interval.start <=
              (*--(*c_previous).end())->read_interval.start))) // i.e. not least one hit outside overlap
            // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, but could also impose no more than 2k hits in overlap
        {
            if (c_previous->size() >= c_current->size()) {
                clusters_of_hits.erase(c_current);
                //cout << "erase current" << endl;
                c_current = c_previous;
            } else {
                clusters_of_hits.erase(c_previous);
                //cout << "erase previous" << endl;
            }
        }
        c_previous = c_current;
    }
    cout << now() << "Now have " << clusters_of_hits.size() << " clusters of hits " << endl;
}

void filter_clusters2(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits, const uint &genome_size) {
    // Sort clusters by size, and filter out those small clusters which are entirely contained in bigger clusters on reads
    cout << now() << "Filter2 the " << clusters_of_hits.size() << " clusters of hits " << endl;
    if (clusters_of_hits.empty()) { return; }

    set<set<MinimizerHitPtr, pComp>, clusterComp_size> clusters_by_size(clusters_of_hits.begin(),
                                                                       clusters_of_hits.end());

    auto it = clusters_by_size.begin();
    std::vector<int> read_v(genome_size, 0);
    //cout << "fill from " << (*(it->begin()))->read_interval.start << " to " << (*--(it->end()))->read_interval.start << endl;
    fill(read_v.begin() + (*(it->begin()))->read_interval.start, read_v.begin() + (*--(it->end()))->read_interval.start,
         1);
    bool contained;
    for (auto it_next = ++clusters_by_size.begin();
         it_next != clusters_by_size.end(); ++it_next) {
        //cout << "read id " << (*(it_next->begin()))->prg_id << endl;
        if ((*(it_next->begin()))->read_id == (*(it->begin()))->read_id) {
            //check if have any 0s in interval of read_v between first and last
            contained = true;
            for (uint i = (*(it_next->begin()))->read_interval.start;
                 i < (*--(it_next->end()))->read_interval.start; ++i) {
                //cout << i << ":" << read_v[i] << "\t";
                if (read_v[i] == 0) {
                    contained = false;
                    //cout << "found unique element at read position " << i << endl;
                    //cout << "fill from " << i << " to " << (*--(it_next->end()))->read_interval.start << endl;
                    fill(read_v.begin() + i, read_v.begin() + (*--(it_next->end()))->read_interval.start, 1);
                    break;
                }

            }
            //cout << endl;
            if (contained) {
                //cout << "erase cluster so clusters_of_hits has size decrease from " << clusters_of_hits.size();
                clusters_of_hits.erase(*it_next);
                //cout << " to " << clusters_of_hits.size() << endl;
            }
        } else {
            //cout << "consider new read" << endl;
            fill(read_v.begin(), read_v.end(), 0);
        }
        ++it;
    }
    cout << now() << "Now have " << clusters_of_hits.size() << " clusters of hits " << endl;
}

void infer_localPRG_order_for_reads(const vector<LocalPRG *> &prgs, MinimizerHits *minimizer_hits, PanGraph *pangraph,
                                    const int max_diff, const uint &genome_size, const float& scale_cluster_size,
                                    const uint min_cluster_size, const uint short_read_length) {
    // this step infers the gene order for a read and adds this to the pangraph
    // by defining clusters of hits, keeping those which are not noise and
    // then adding the inferred gene ordering

    minimizer_hits->sort();
    if (minimizer_hits->hits.empty()) { return; }

    set<set<MinimizerHitPtr, pComp>, clusterComp> clusters_of_hits;
    define_clusters(clusters_of_hits, prgs, minimizer_hits, max_diff, scale_cluster_size,
                    min_cluster_size, short_read_length);

    filter_clusters(clusters_of_hits);
    //filter_clusters2(clusters_of_hits, genome_size);

    // Add inferred order to pangraph    
    if (clusters_of_hits.empty()) { return; }
    // to do this consider pairs of clusters in turn
    auto c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, prgs[(*(*c_previous).begin())->prg_id]->name,
                       (*(*c_previous).begin())->read_id, *c_previous);
    //cout << "nodes on read " << (*(*c_previous).begin())->read_id << " : "
    //     << prgs[(*(*c_previous).begin())->prg_id]->name;
    for (auto c_current = ++clusters_of_hits.begin();
         c_current != clusters_of_hits.end(); ++c_current) {
        if ((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) {
            pangraph->add_node((*(*c_current).begin())->prg_id, prgs[(*(*c_current).begin())->prg_id]->name,
                               (*(*c_current).begin())->read_id, *c_current);
            if ((*(*c_previous).begin())->prg_id != (*(*c_current).begin())->prg_id or
                (*(*c_previous).begin())->strand != (*(*c_current).begin())->strand) {
                pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id,
                                   (*(*c_previous).begin())->strand + 2 * (*(*c_current).begin())->strand,
                                   (*(*c_current).begin())->read_id);
                //cout << "\t" << prgs[(*(*c_current).begin())->prg_id]->name;
            }
            c_previous = c_current;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id) {
            // if we just started looking at hits for a new read, add the first cluster
            //cout << endl << "nodes on read " << (*(*c_current).begin())->read_id << " : "
            //     << prgs[(*(*c_current).begin())->prg_id]->name;
            pangraph->add_node((*(*c_current).begin())->prg_id, prgs[(*(*c_current).begin())->prg_id]->name,
                               (*(*c_current).begin())->read_id, *c_current);
            c_previous = c_current;
        }
    }
    cout << endl;
}

void pangraph_from_read_file(const string &filepath, MinimizerHits *mh, PanGraph *pangraph, Index *idx,
                             const vector<LocalPRG *> &prgs, const uint32_t w, const uint32_t k,
                             const int max_diff, const float& e_rate,
                             const uint min_cluster_size, const uint genome_size,
                             const bool illumina)
{
    string name, read, line;
    uint64_t covg = 0;
    float scale_cluster_size = 0.75 / exp(e_rate * k);
    uint short_read_length = std::numeric_limits<uint>::max();
    Seq *s;
    s = new Seq(0, "null", "", w, k);
    if (s == nullptr) {
        cerr << "Failed to create new Seq, something must be dying " << endl;
        exit(EXIT_FAILURE);
    }
    uint32_t id = 0;

    ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line.empty() || line[0] == '>' || line[0] == '@') {
                if (!read.empty()) // ok we'll allow reads with no name, removed
                {
                    s->initialize(id, name, read, w, k);
                    covg += s->seq.length();
                    if (illumina == true and short_read_length == std::numeric_limits<uint>::max())
                    {
                        short_read_length = s->seq.length();
                    }
                    //cout << now() << "Add read hits" << endl;
                    add_read_hits(s, mh, idx);
                    id++;
                }
                name.clear();
                read.clear();
                if (!line.empty()) {
                    name = line.substr(1);
                }
            } else if (line[0] == '+') {
                //skip this line and the qual score line
                getline(myfile, line);
            } else {
                read += line;
            }
        }
        // and last entry
        if (!read.empty()) // allow reads with no name
        {
            //cout << now() << "Found read " << name << endl;
            s->initialize(id, name, read, w, k);
            covg += s->seq.length();
            if (illumina == true and short_read_length == std::numeric_limits<uint>::max())
            {
                short_read_length = s->seq.length();
            }
            //cout << now() << "Add read hits" << endl;
            add_read_hits(s, mh, idx);
        }
        covg = covg / genome_size;
        cout << now() << "Estimated coverage: " << covg << endl;
        //cout << "Number of reads found: " << id+1 << endl;
        cout << now() << "Infer gene orders and add to PanGraph" << endl;
        infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, genome_size, scale_cluster_size,
                                       min_cluster_size, short_read_length);
        cout << now() << "Pangraph has " << pangraph->nodes.size() << " nodes" << endl;
        pangraph->clean(covg);
        cout << now() << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes" << endl;
        delete s;
        myfile.close();
    } else {
        cerr << "Unable to open read file " << filepath << endl;
        exit(EXIT_FAILURE);
    }
}

void update_localPRGs_with_hits(PanGraph *pangraph,
                                const vector<LocalPRG *> &prgs) //, const uint32_t k, const float& e_rate, bool output_p_dist)
{
    pangraph->add_hits_to_kmergraphs(prgs);
}

