#include <cstring>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cctype>
#include <set>
#include <memory>
#include <ctime>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "utils.h"
#include "seq.h"
#include "localPRG.h"
#include "pangenome/pangraph.h"
#include "pangenome/panread.h"
#include "noise_filtering.h"
#include "minihit.h"
#include "fastaq_handler.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

string now() {
    time_t now;
    string dt;

    now = time(nullptr);
    dt = ctime(&now);
    return dt.substr(0, dt.length() - 1) + " ";
}

void make_dir(const string &dirpath) {
    if (dirpath == ".")
        return;

    const boost::filesystem::path dir(dirpath);
    if (boost::filesystem::exists(dir)) {
        return;
    }

    assert(boost::filesystem::create_directories(dir) || assert_msg("Failed to make directories " << dir));
}

string int_to_string(const int number) {
    stringstream ss;
    ss << setw(2) << setfill('0') << number;
    return ss.str();
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

    for (uint32_t m = n; m != n - k1 - k2; --m) {
        total += log(m);
    }

    for (uint32_t m = 1; m < k1; ++m) {
        total -= log(m + 1);
    }

    for (uint32_t m = 1; m < k2; ++m) {
        total -= log(m + 1);
    }

    return total;
}

void read_prg_file(std::vector<std::shared_ptr<LocalPRG>> &prgs,
                   const std::string &filepath) {
    BOOST_LOG_TRIVIAL(debug) << "Loading PRGs from file " << filepath;

    uint32_t id = 0;

    FastaqHandler fh(filepath);
    while (!fh.eof()) {
        fh.get_next();
        if (fh.name.empty() or fh.read.empty())
            continue;
        auto s = std::make_shared<LocalPRG>(LocalPRG(id, fh.name, fh.read));
        if (s != nullptr) {
            prgs.push_back(s);
            id++;
        } else {
            cerr << "Failed to make LocalPRG for " << fh.name << endl;
            exit(1);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Number of LocalPRGs read: " << prgs.size();
}

void load_PRG_kmergraphs(std::vector<std::shared_ptr<LocalPRG>> &prgs, const uint32_t &w, const uint32_t &k,
                         const string &prgfile) {
    BOOST_LOG_TRIVIAL(debug) << "Loading kmer_prgs from files";
    string prefix = "";
    size_t pos = prgfile.find_last_of("/");
    if (pos != std::string::npos) {
        prefix += prgfile.substr(0, pos);
        prefix += "/";
    }
    //cout << "prefix for kmerprgs dir is " << prefix << endl;

    auto dir_num = 0;
    string dir;
    for (const auto &prg: prgs) {
        //cout << "Load kmergraph for " << prg->name << endl;
        if (prg->id % 4000 == 0) {
            dir = prefix + "kmer_prgs/" + int_to_string(dir_num + 1);
            dir_num++;
            boost::filesystem::path p(dir);
            if (not boost::filesystem::exists(p))
                dir = prefix + "kmer_prgs";
        }
        prg->kmer_prg.load(dir + "/" + prg->name + ".k" + to_string(k) + ".w" + to_string(w) + ".gfa");
    }
}

void load_vcf_refs_file(const string &filepath, VCFRefs &vcf_refs) {
    BOOST_LOG_TRIVIAL(debug) << "Loading VCF refs from file " << filepath;

    FastaqHandler fh(filepath);
    while (!fh.eof()) {
        fh.get_next();
        if (!fh.name.empty() && !fh.read.empty()) {
            vcf_refs[fh.name] = fh.read;
        }
    }
}

//void add_read_hits(const uint32_t id, const string& name, const string& seq, MinimizerHits* hits, Index* idx, const uint32_t w, const uint32_t k)
void add_read_hits(Seq *s, MinimizerHits *hits, Index *idx) {
    //cout << now() << "Search for hits for read " << s->name << " which has sketch size " << s->sketch.size() << " against index of size " << idx->minhash.size() << endl;
    uint32_t hit_count = 0;
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and adds hits to a global MinimizerHits object
    //Seq s(id, name, seq, w, k);
    for (auto it = s->sketch.begin(); it != s->sketch.end(); ++it) {
        if (idx->minhash.find((*it).kmer) != idx->minhash.end()) {
            for (uint32_t j = 0; j != idx->minhash[(*it).kmer]->size(); ++j) {
                hits->add_hit(s->id, *it, &(idx->minhash[(*it).kmer]->operator[](j)));
                hit_count += 1;
            }
            //} else {
            //    cout << "did not find minimizer " << (*it).kmer << " in index" << endl;
        }
    }
    //hits->sort();
    //cout << now() << "Found " << hit_count << " hits found for read " << s->name << " so size of MinimizerHits is now "
    //     << hits->hits.size() + hits->uhits.size() << endl;
}

void define_clusters(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits,
                     const std::vector<std::shared_ptr<LocalPRG>> &prgs,
                     MinimizerHits *minimizer_hits, const int max_diff,
                     const float &fraction_kmers_required_for_cluster,
                     const uint32_t min_cluster_size, const uint32_t expected_number_kmers_in_short_read_sketch) {
    BOOST_LOG_TRIVIAL(debug) << "Define clusters of hits from the " << minimizer_hits->hits.size() << " hits";

    if (minimizer_hits->hits.empty()) { return; }

    // A cluster of hits should match same localPRG, each hit not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    auto mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHitPtr, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    uint32_t length_based_threshold;
    for (auto mh_current = ++minimizer_hits->hits.begin();
         mh_current != minimizer_hits->hits.end(); ++mh_current) {
        if ((*mh_current)->read_id != (*mh_previous)->read_id or
            (*mh_current)->prg_id != (*mh_previous)->prg_id or
            (*mh_current)->strand != (*mh_previous)->strand or
            (abs((int) (*mh_current)->read_start_position - (int) (*mh_previous)->read_start_position)) > max_diff) {
            // keep clusters which cover at least 3/4 the expected number of minihits
            length_based_threshold = min(prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length(),
                                         expected_number_kmers_in_short_read_sketch) *
                                     fraction_kmers_required_for_cluster;
            /*cout << "gene length " << prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()
                 << " short read length: " << expected_number_kmers_in_short_read_sketch
                 << " scale cluster size: " << fraction_kmers_required_for_cluster
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
                                 expected_number_kmers_in_short_read_sketch) * fraction_kmers_required_for_cluster;
    /*cout << "gene length " << prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()
         << " short read length: " << expected_number_kmers_in_short_read_sketch
         << " scale cluster size: " << fraction_kmers_required_for_cluster
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

    BOOST_LOG_TRIVIAL(debug) << "Found " << clusters_of_hits.size() << " clusters of hits";
}

void filter_clusters(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits) {
    // Next order clusters, choose between those that overlap by too much
    BOOST_LOG_TRIVIAL(debug) << "Filter the " << clusters_of_hits.size() << " clusters of hits";
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
             ((*--(*c_current).end())->read_start_position <=
              (*--(*c_previous).end())->read_start_position))) // i.e. not least one hit outside overlap
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
    BOOST_LOG_TRIVIAL(debug) << "Now have " << clusters_of_hits.size() << " clusters of hits";
}

void filter_clusters2(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits, const uint32_t &genome_size) {
    // Sort clusters by size, and filter out those small clusters which are entirely contained in bigger clusters on reads
    BOOST_LOG_TRIVIAL(debug) << "Filter2 the " << clusters_of_hits.size() << " clusters of hits";
    if (clusters_of_hits.empty()) { return; }

    set<set<MinimizerHitPtr, pComp>, clusterComp_size> clusters_by_size(clusters_of_hits.begin(),
                                                                        clusters_of_hits.end());

    auto it = clusters_by_size.begin();
    std::vector<int> read_v(genome_size, 0);
    //cout << "fill from " << (*(it->begin()))->read_start_position << " to " << (*--(it->end()))->read_start_position << endl;
    fill(read_v.begin() + (*(it->begin()))->read_start_position, read_v.begin() + (*--(it->end()))->read_start_position,
         1);
    bool contained;
    for (auto it_next = ++clusters_by_size.begin();
         it_next != clusters_by_size.end(); ++it_next) {
        //cout << "read id " << (*(it_next->begin()))->prg_id << endl;
        if ((*(it_next->begin()))->read_id == (*(it->begin()))->read_id) {
            //check if have any 0s in interval of read_v between first and last
            contained = true;
            for (uint32_t i = (*(it_next->begin()))->read_start_position;
                 i < (*--(it_next->end()))->read_start_position; ++i) {
                //cout << i << ":" << read_v[i] << "\t";
                if (read_v[i] == 0) {
                    contained = false;
                    //cout << "found unique element at read position " << i << endl;
                    //cout << "fill from " << i << " to " << (*--(it_next->end()))->read_start_position << endl;
                    fill(read_v.begin() + i, read_v.begin() + (*--(it_next->end()))->read_start_position, 1);
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
    BOOST_LOG_TRIVIAL(debug) << "Now have " << clusters_of_hits.size() << " clusters of hits";
}

void
add_clusters_to_pangraph(set<set<MinimizerHitPtr, pComp>, clusterComp> &clusters_of_hits, pangenome::Graph *pangraph,
                         const std::vector<std::shared_ptr<LocalPRG>> &prgs) {
    BOOST_LOG_TRIVIAL(debug) << "Add inferred order to PanGraph";
    if (clusters_of_hits.empty()) { return; }

    // to do this consider pairs of clusters in turn
    for (auto cluster: clusters_of_hits) {
        pangraph->add_node((*cluster.begin())->prg_id,
                           prgs[(*cluster.begin())->prg_id]->name,
                           (*cluster.begin())->read_id,
                           cluster);
    }
}

void infer_localPRG_order_for_reads(const std::vector<std::shared_ptr<LocalPRG>> &prgs, MinimizerHits *minimizer_hits,
                                    pangenome::Graph *pangraph,
                                    const int max_diff, const uint32_t &genome_size,
                                    const float &fraction_kmers_required_for_cluster,
                                    const uint32_t min_cluster_size,
                                    const uint32_t expected_number_kmers_in_short_read_sketch) {
    // this step infers the gene order for a read and adds this to the pangraph
    // by defining clusters of hits, keeping those which are not noise and
    // then adding the inferred gene ordering

    minimizer_hits->sort();
    if (minimizer_hits->hits.empty()) { return; }

    set<set<MinimizerHitPtr, pComp>, clusterComp> clusters_of_hits;
    define_clusters(clusters_of_hits, prgs, minimizer_hits, max_diff, fraction_kmers_required_for_cluster,
                    min_cluster_size, expected_number_kmers_in_short_read_sketch);

    minimizer_hits->clear();

    filter_clusters(clusters_of_hits);
    //filter_clusters2(clusters_of_hits, genome_size);
    add_clusters_to_pangraph(clusters_of_hits, pangraph, prgs);
}

uint32_t pangraph_from_read_file(const string &filepath,
                                 MinimizerHits *mh,
                                 pangenome::Graph *pangraph,
                                 Index *idx,
                                 const std::vector<std::shared_ptr<LocalPRG>> &prgs,
                                 const uint32_t w,
                                 const uint32_t k,
                                 const int max_diff,
                                 const float &e_rate,
                                 const uint32_t min_cluster_size,
                                 const uint32_t genome_size,
                                 const bool illumina,
                                 const bool clean,
                                 const uint32_t max_covg) {
    uint64_t covg = 0;
    float fraction_kmers_required_for_cluster = 0.75 / exp(e_rate * k);
    uint32_t expected_number_kmers_in_short_read_sketch = std::numeric_limits<uint32_t>::max();
    Seq *s;
    s = new Seq(0, "null", "", w, k);
    if (s == nullptr) {
        cerr << "Failed to create new Seq, something must be dying " << endl;
        exit(EXIT_FAILURE);
    }
    uint32_t id = 0;

    FastaqHandler fh(filepath);
    while (!fh.eof()) {
        fh.get_next();
        s->initialize(id, fh.name, fh.read, w, k);
        if (!s->sketch.empty()) {
            covg += s->seq.length();
            if (covg / genome_size > max_covg) {
                BOOST_LOG_TRIVIAL(warning) << "Stop reading readfile as have reached max coverage";
                break;
            }
        } else {
            id++;
            continue;
        }
        if (illumina and expected_number_kmers_in_short_read_sketch == std::numeric_limits<uint32_t>::max()) {
            assert(w != 0);
            expected_number_kmers_in_short_read_sketch = s->seq.length() * 2 / w;
        }
        //cout << now() << "Add read hits" << endl;
        add_read_hits(s, mh, idx);
        id++;
        if (id > 10000000) {
            BOOST_LOG_TRIVIAL(debug) << "Stop reading readfile as have reached 10,000,000 reads";
            break;
        }

        if (mh->uhits.size() > 90000) {
            BOOST_LOG_TRIVIAL(debug) << "Infer gene orders and add to pangenome::Graph";
            pangraph->reserve_num_reads(id);
            infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, genome_size,
                                           fraction_kmers_required_for_cluster, min_cluster_size,
                                           expected_number_kmers_in_short_read_sketch);
        }
    }
    delete s;
    BOOST_LOG_TRIVIAL(debug) << "Found " << id << " reads";

    BOOST_LOG_TRIVIAL(debug) << "Infer gene orders and add to pangenome::Graph";
    pangraph->reserve_num_reads(id);
    infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, genome_size, fraction_kmers_required_for_cluster,
                                   min_cluster_size, expected_number_kmers_in_short_read_sketch);

    BOOST_LOG_TRIVIAL(debug) << "Pangraph has " << pangraph->nodes.size() << " nodes";

    covg = covg / genome_size;
    BOOST_LOG_TRIVIAL(debug) << "Estimated coverage: " << covg;

    if (illumina and clean) {
        clean_pangraph_with_debruijn_graph(pangraph, 2, 1, illumina);
        BOOST_LOG_TRIVIAL(debug) << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes";
    } else if (clean) {
        clean_pangraph_with_debruijn_graph(pangraph, 3, 1, illumina);
        BOOST_LOG_TRIVIAL(debug) << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes";
    }

    return covg;
}

void update_localPRGs_with_hits(pangenome::Graph *pangraph,
                                const std::vector<std::shared_ptr<LocalPRG>> &prgs) //, const uint32_t k, const float& e_rate, bool output_p_dist)
{
    pangraph->add_hits_to_kmergraphs(prgs);
}

