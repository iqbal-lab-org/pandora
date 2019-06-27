#include <sstream>
#include <fstream>
#include <cassert>
#include <limits>
#include <cstdio>      /* NULL */
#include <cstdlib>     /* srand, rand */
#include <cmath>
#include <cstdint>

#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/log/trivial.hpp>

#include "utils.h"
#include "kmernode.h"
#include "kmergraphwithcoverage.h"
#include "localPRG.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


using namespace prg;

void KmerGraphWithCoverage::set_exp_depth_covg(const uint32_t edp) {
    assert(edp > 0);
    exp_depth_covg = edp;
}

void KmerGraphWithCoverage::set_binomial_parameter_p(const float e_rate) {
    BOOST_LOG_TRIVIAL(debug) << "Set p in kmergraph";
    assert(kmer_prg->k != 0);
    assert(0 < e_rate and e_rate < 1);
    binomial_parameter_p = 1 / exp(e_rate * kmer_prg->k);
    //cout << "using p: " << p << endl;
}

void KmerGraphWithCoverage::increment_covg(uint32_t node_id, bool strand, uint32_t sample_id) {
    assert(this->nodeIndex2SampleCoverage[node_id].size() > sample_id);

    //get a pointer to the value we want to increment
    uint16_t* coverage_ptr = nullptr;
    if (strand)
        coverage_ptr = &(this->nodeIndex2SampleCoverage[node_id][sample_id].first);
    else
        coverage_ptr = &(this->nodeIndex2SampleCoverage[node_id][sample_id].second);

    if ((*coverage_ptr) < UINT16_MAX) //checks if it is safe to increase the coverage (no overflow)
        ++(*coverage_ptr);
}

uint32_t KmerGraphWithCoverage::get_covg(uint32_t node_id, bool strand, uint32_t sample_id) const {

    if (this->nodeIndex2SampleCoverage[node_id].size() <= sample_id)
        return 0;

    if (strand)
        return (uint32_t)(this->nodeIndex2SampleCoverage[node_id][sample_id].first);
    else
        return (uint32_t)(this->nodeIndex2SampleCoverage[node_id][sample_id].second);
}

void KmerGraphWithCoverage::set_covg(uint32_t node_id, uint16_t value, bool strand, uint32_t sample_id) {
    assert(this->nodeIndex2SampleCoverage[node_id].size() > sample_id);
    if (strand)
        this->nodeIndex2SampleCoverage[node_id][sample_id].first = value;
    else
        this->nodeIndex2SampleCoverage[node_id][sample_id].second = value;
}



void KmerGraphWithCoverage::set_negative_binomial_parameters(const float &nbin_prob, const float &nb_fail) {
    if (nbin_prob == 0 and nb_fail == 0)
        return;
    //qcout << "set nb" << endl;
    assert((negative_binomial_parameter_p > 0 and negative_binomial_parameter_p < 1) || assert_msg("nb_p " << negative_binomial_parameter_p << " was not set in kmergraph"));
    assert(negative_binomial_parameter_r > 0 || assert_msg("nb_r was not set in kmergraph"));
    negative_binomial_parameter_p += nbin_prob;
    negative_binomial_parameter_r += nb_fail;
}


float KmerGraphWithCoverage::nbin_prob(uint32_t node_id, const uint32_t &sample_id) {
    auto k = get_covg(node_id, 0, sample_id) + get_covg(node_id, 1, sample_id);
    float return_prob = log(pdf(boost::math::negative_binomial(negative_binomial_parameter_r, negative_binomial_parameter_p), k));
    return_prob = std::max(return_prob, std::numeric_limits<float>::lowest() / 1000);
    return return_prob;
}

float KmerGraphWithCoverage::lin_prob(uint32_t node_id, const uint32_t &sample_id) {
    assert(num_reads != 0);
    auto k = get_covg(node_id, 0, sample_id) + get_covg(node_id, 1, sample_id);
    return log(float(k)/num_reads);
}

float KmerGraphWithCoverage::bin_prob(uint32_t node_id, const uint32_t &sample_id) {
    assert(num_reads != 0);
    return bin_prob(node_id, num_reads, sample_id);
}

float KmerGraphWithCoverage::bin_prob(const uint32_t &node_id, const uint32_t &num, const uint32_t &sample_id) {
    assert(binomial_parameter_p != 1);
    assert(node_id < kmer_prg->nodes.size());
    #ifndef NDEBUG
        //TODO: check if tests must be updated or not due to this (I think not - sorted_nodes is always sorted)
        //if this is added, some tests bug, since it was not executed before...
        //check();
    #endif

    uint32_t sum_coverages = get_covg(node_id, 0, sample_id) + get_covg(node_id, 1, sample_id);

    float return_prob;
    if (node_id == (*(kmer_prg->sorted_nodes.begin()))->id or node_id == (*(kmer_prg->sorted_nodes.rbegin()))->id) {
        return_prob = 0; // is really undefined
    } else if (sum_coverages > num) {
        // under model assumptions this can't happen, but it inevitably will, so bodge
        return_prob = lognchoosek2(sum_coverages,
                           get_covg(node_id, 0, sample_id),
                           get_covg(node_id, 1, sample_id))
              + sum_coverages * log(binomial_parameter_p / 2);
        // note this may give disadvantage to repeat kmers
    } else {
        return_prob = lognchoosek2(num,
                           get_covg(node_id, 0, sample_id),
                           get_covg(node_id, 1, sample_id))
              + sum_coverages * log(binomial_parameter_p / 2)
              + (num - sum_coverages) * log(1 - binomial_parameter_p);
    }
    return return_prob;
}

float KmerGraphWithCoverage::get_prob(const std::string& prob_model, const uint32_t &node_id, const uint32_t &sample_id){
    if (prob_model == "nbin")
    {
        return nbin_prob(node_id, sample_id);
    } else if (prob_model == "bin") {
        assert(binomial_parameter_p < 1 || assert_msg("binomial_parameter_p was not set in kmergraph"));
        assert(num_reads > 0 || assert_msg("num_reads was not set in kmergraph"));
        return bin_prob(node_id, sample_id);
    } else if (prob_model == "lin") {
        return lin_prob(node_id, sample_id);
    } else {
        BOOST_LOG_TRIVIAL(warning) << "Invalid probability model for kmer coverage distribution: should be nbin, bin or lin";
        exit(1);
    }
}

bool KmerGraphWithCoverage::coverage_is_zeroes(const uint32_t& sample_id){
    bool all_zero = true;
    for (const auto &node_ptr : kmer_prg->nodes) {
        if (get_covg(node_ptr->id, 0, sample_id) + get_covg(node_ptr->id, 1, sample_id) > 0) {
            BOOST_LOG_TRIVIAL(debug) << "Found non-zero coverage in kmer graph";
            all_zero = false;
            break;
        }
    }
    if (all_zero) {
        BOOST_LOG_TRIVIAL(debug) << "ALL ZEROES in kmer graph coverages";
    }
    return all_zero;
}

float KmerGraphWithCoverage::find_max_path(std::vector<KmerNodePtr> &maxpath,
                               const std::string& prob_model,
                               const uint32_t &max_num_kmers_to_average,
                               const uint32_t &sample_id) {
    //TODO: FIX THIS INNEFICIENCY I INTRODUCED
    const std::vector<KmerNodePtr> sorted_nodes(this->kmer_prg->sorted_nodes.begin(),
                                                this->kmer_prg->sorted_nodes.end());
    this->kmer_prg->check();

    // also check not all 0 covgs
    auto coverages_all_zero = coverage_is_zeroes(sample_id);
    if (coverages_all_zero)
        return std::numeric_limits<float>::lowest();

    // create vectors to hold the intermediate values
    std::vector<float> max_sum_of_log_probs_from_node(sorted_nodes.size(), 0);
    std::vector<uint32_t> length_of_maxpath_from_node(sorted_nodes.size(), 0);
    std::vector<uint32_t> prev_node_along_maxpath(sorted_nodes.size(), sorted_nodes.size() - 1);
    float max_mean;
    int max_length;
    const float tolerance = 0.000001;

    for (uint32_t j = sorted_nodes.size() - 1; j != 0; --j) {
        max_mean = std::numeric_limits<float>::lowest();
        max_length = 0; // tie break with longest kmer path
        const auto & current_node = sorted_nodes[j - 1];
        for (uint32_t i = 0; i != current_node->out_nodes.size(); ++i) {
            const auto & considered_outnode = current_node->out_nodes[i].lock();
            bool is_terminus_and_most_likely = considered_outnode->id == sorted_nodes.back()->id
                                               and thresh > max_mean + tolerance;
            bool avg_log_likelihood_is_most_likely = max_sum_of_log_probs_from_node[considered_outnode->id] /
                                                     length_of_maxpath_from_node[considered_outnode->id]
                                                     > max_mean + tolerance;
            bool avg_log_likelihood_is_close_to_most_likely = max_mean
                                                              - max_sum_of_log_probs_from_node[considered_outnode->id] /
                                                                length_of_maxpath_from_node[considered_outnode->id]
                                                              <= tolerance;
            bool is_longer_path = length_of_maxpath_from_node[considered_outnode->id] > max_length;

            if (is_terminus_and_most_likely or avg_log_likelihood_is_most_likely
                or ( avg_log_likelihood_is_close_to_most_likely and is_longer_path)) {
                max_sum_of_log_probs_from_node[current_node->id] = get_prob(prob_model,current_node->id, sample_id)
                                      + max_sum_of_log_probs_from_node[considered_outnode->id];
                length_of_maxpath_from_node[current_node->id] = 1 + length_of_maxpath_from_node[considered_outnode->id];
                prev_node_along_maxpath[current_node->id] = considered_outnode->id;

                if (length_of_maxpath_from_node[current_node->id] > max_num_kmers_to_average){
                    uint32_t prev_node = prev_node_along_maxpath[current_node->id];
                    for (uint step=0; step<max_num_kmers_to_average; step++){
                        prev_node = prev_node_along_maxpath[prev_node];
                    }
                    max_sum_of_log_probs_from_node[current_node->id] -= get_prob(prob_model,sorted_nodes[prev_node]->id,
                                                                                 sample_id);
                    length_of_maxpath_from_node[current_node->id] -= 1;
                    assert(length_of_maxpath_from_node[current_node->id] == max_num_kmers_to_average);
                }

                if (considered_outnode->id != sorted_nodes.back()->id) {
                    max_mean = max_sum_of_log_probs_from_node[considered_outnode->id] /
                            length_of_maxpath_from_node[considered_outnode->id];
                    max_length = length_of_maxpath_from_node[considered_outnode->id];
                } else {
                    max_mean = thresh;
                }
            }
        }
    }

    // extract path
    uint32_t prev_node = prev_node_along_maxpath[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1) {
        maxpath.push_back(this->kmer_prg->nodes[prev_node]);
        prev_node = prev_node_along_maxpath[prev_node];

        if (maxpath.size() > 1000000){
            BOOST_LOG_TRIVIAL(warning) << "I think I've found an infinite loop - is something wrong with this kmergraph?";
            exit(1);
        }
    }

    assert(length_of_maxpath_from_node[0] > 0 || assert_msg("found no path through kmer prg"));
    return prob_path(maxpath, sample_id, prob_model);
}


/*
std::vector<std::vector<KmerNodePtr>> KmerGraph::find_max_paths(uint32_t num,
                                                                const uint32_t &sample_id) {

    // save original coverges so can put back at the end
    std::vector<uint32_t> original_covgs0, original_covgs1;
    for (uint32_t i = 0; i != nodes.size(); ++i) {
        original_covgs0.push_back(nodes[i]->get_covg(0, sample_id));
        original_covgs1.push_back(nodes[i]->get_covg(1, sample_id));
    }

    // find num max paths
    //cout << "expected covg " << (uint)(p*num_reads/num) << endl;
    std::vector<std::vector<KmerNodePtr>> paths;
    std::vector<KmerNodePtr> maxpath;
    find_max_path(maxpath, <#initializer#>);
    //uint min_covg;
    paths.push_back(maxpath);

    while (paths.size() < num) {
        for (uint32_t i = 0; i != maxpath.size(); ++i) {
            maxpath[i]->covg[0] -= std::min(maxpath[i]->covg[0], (uint32_t) (p * num_reads / num));
            maxpath[i]->covg[1] -= std::min(maxpath[i]->covg[1], (uint32_t) (p * num_reads / num));
        }
        maxpath.clear();
        find_max_path(maxpath, <#initializer#>);
        paths.push_back(maxpath);
    }

    // put covgs back
    for (uint32_t i = 0; i != nodes.size(); ++i) {
        nodes[i]->covg[0] = original_covgs0[i];
        nodes[i]->covg[1] = original_covgs1[i];
    }

    return paths;
}
 */

std::vector<std::vector<KmerNodePtr>> KmerGraphWithCoverage::get_random_paths(uint32_t num_paths) {
    // find a random path through kmergraph picking ~uniformly from the outnodes at each point
    std::vector<std::vector<KmerNodePtr>> rpaths;
    std::vector<KmerNodePtr> rpath;
    uint32_t i;

    time_t now;
    now = time(nullptr);
    srand((unsigned int) now);

    if (!kmer_prg->nodes.empty()) {
        for (uint32_t j = 0; j != num_paths; ++j) {
            i = rand() % kmer_prg->nodes[0]->out_nodes.size();
            rpath.push_back(kmer_prg->nodes[0]->out_nodes[i].lock());
            while (rpath.back() != kmer_prg->nodes[kmer_prg->nodes.size() - 1]) {
                if (rpath.back()->out_nodes.size() == 1) {
                    rpath.push_back(rpath.back()->out_nodes[0].lock());
                } else {
                    i = rand() % rpath.back()->out_nodes.size();
                    rpath.push_back(rpath.back()->out_nodes[i].lock());
                }
            }
            rpath.pop_back();
            rpaths.push_back(rpath);
            rpath.clear();
        }
    }
    return rpaths;
}

float KmerGraphWithCoverage::prob_path(const std::vector<KmerNodePtr> &kpath,
                           const uint32_t &sample_id, const std::string& prob_model) {
    float return_prob_path = 0;
    for (uint32_t i = 0; i != kpath.size(); ++i) {
        return_prob_path += get_prob(prob_model, kpath[i]->id, sample_id);
    }
    uint32_t len = kpath.size();
    if (kpath[0]->path.length() == 0) {
        len -= 1;
    }
    if (kpath.back()->path.length() == 0) {
        len -= 1;
    }
    if (len == 0) {
        len = 1;
    }
    return return_prob_path / len;
}

/*
float KmerGraph::prob_paths(const std::vector<std::vector<KmerNodePtr>> &kpaths) {
    if (kpaths.empty()) {
        return 0; // is this the correct default?
    }

    // collect how much coverage we expect on each node from these paths
    std::vector<uint32_t> path_node_covg(nodes.size(), 0);
    for (uint32_t i = 0; i != kpaths.size(); ++i) {
        for (uint32_t j = 0; j != kpaths[i].size(); ++j) {
            path_node_covg[kpaths[i][j]->id] += 1;
        }
    }

    // now calculate max likelihood assuming independent paths
    float ret_p = 0;
    uint32_t len = 0;
    for (uint32_t i = 0; i != path_node_covg.size(); ++i) {
        if (path_node_covg[i] > 0) {
            //cout << "prob of node " << nodes[i]->id << " which has path covg " << path_node_covg[i] << " and so we expect to see " << num_reads*path_node_covg[i]/kpaths.size() << " times IS " << prob(nodes[i]->id, num_reads*path_node_covg[i]/kpaths.size()) << endl;
                ret_p += prob(nodes[i]->id, num_reads * path_node_covg[i] / kpaths.size(), <#initializer#>);
            if (nodes[i]->path.length() > 0) {
                len += 1;
            }
        }
    }

    if (len == 0) {
        len = 1;
    }

    //cout << "len " << len << endl;
    return ret_p / len;
}
*/

void KmerGraphWithCoverage::save_covg_dist(const std::string &filepath) {
    std::ofstream handle;
    handle.open(filepath);

    for (const auto &kmer_node_ptr: kmer_prg->nodes) {
        const KmerNode &kmer_node = *kmer_node_ptr;

        uint32_t sample_id = 0;
        for (const auto &sample_coverage: nodeIndex2SampleCoverage[kmer_node.id]) {
            handle << kmer_node.id << " "
                   << sample_id << " "
                   << sample_coverage.first << " "
                   << sample_coverage.second;

            sample_id++;
        }
    }
    handle.close();
}

//save the KmerGraph as gfa
//TODO: THIS SHOULD BE RECODED, WE ARE DUPLICATING CODE HERE (SEE KmerGraph::save())!!!
void KmerGraphWithCoverage::save(const std::string &filepath, const std::shared_ptr<LocalPRG> localprg) {
    uint32_t sample_id = 0;

    std::ofstream handle;
    handle.open(filepath);
    if (handle.is_open()) {
        handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
        for (const auto &c : kmer_prg->nodes) {
            handle << "S\t" << c->id << "\t";

            if (localprg != nullptr) {
                handle << localprg->string_along_path(c->path);
            } else {
                handle << c->path;
            }

            handle << "\tFC:i:" << get_covg(c->id, 0, sample_id) << "\t" << "\tRC:i:"
                   << get_covg(c->id, 1, sample_id) << std::endl;//"\t" << (unsigned)nodes[i].second->num_AT << endl;

            for (uint32_t j = 0; j < c->out_nodes.size(); ++j) {
                handle << "L\t" << c->id << "\t+\t" << c->out_nodes[j].lock()->id << "\t+\t0M" << std::endl;
            }
        }
        handle.close();
    } else {
        std::cerr << "Unable to open kmergraph file " << filepath << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//TODO: THIS SHOULD BE RECODED, WE ARE DUPLICATING CODE HERE (SEE KmerGraph::load())!!!
//TODO: remove this method?
void KmerGraphWithCoverage::load(const std::string &filepath) {
    //TODO: this might be dangerous, recode this?
    auto kmer_prg = const_cast<KmerGraph*>(this->kmer_prg);
    kmer_prg->clear();
    uint32_t sample_id = 0;

    std::string line;
    std::vector<std::string> split_line;
    std::stringstream ss;
    uint32_t id = 0, from, to;
    uint16_t covg;
    prg::Path p;
    uint32_t num_nodes = 0;

    std::ifstream myfile(filepath);
    if (myfile.is_open()) {

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = std::stoi(split_line[1]);
                num_nodes = std::max(num_nodes, id);
            }
        }
        myfile.clear();
        myfile.seekg(0, myfile.beg);
        kmer_prg->nodes.reserve(num_nodes);
        std::vector<uint16_t> outnode_counts(num_nodes + 1, 0), innode_counts(num_nodes + 1, 0);

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = stoi(split_line[1]);
                ss << split_line[2];
                char c = ss.peek();
                assert (isdigit(c) or assert_msg("Cannot read in this sort of kmergraph GFA as it does not label nodes "
                                                 "with their PRG path"));
                ss >> p;
                ss.clear();
                //add_node(p);
                KmerNodePtr n = std::make_shared<KmerNode>(id, p);
                assert(n != nullptr);
                assert(id == kmer_prg->nodes.size() or num_nodes - id == kmer_prg->nodes.size() or
                       assert_msg("id " << id << " != " << kmer_prg->nodes.size() << " nodes.size() for kmergraph "));
                kmer_prg->nodes.push_back(n);
                kmer_prg->sorted_nodes.insert(n);
                if (kmer_prg->k == 0 and p.length() > 0) {
                    kmer_prg->k = p.length();
                }
                covg = (uint16_t)(stoul(split(split_line[3], "FC:i:")[0]));
                set_covg(n->id, covg, 0, sample_id);
                covg = (uint16_t)(stoul(split(split_line[4], "RC:i:")[0]));
                set_covg(n->id, covg, 1, sample_id);
                if (split_line.size() >= 6) {
                    n->num_AT = std::stoi(split_line[5]);
                }
            } else if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                assert(stoi(split_line[1]) < (int) outnode_counts.size() or
                       assert_msg(stoi(split_line[1]) << ">=" << outnode_counts.size()));
                assert(stoi(split_line[3]) < (int) innode_counts.size() or
                       assert_msg(stoi(split_line[3]) << ">=" << innode_counts.size()));
                outnode_counts[stoi(split_line[1])] += 1;
                innode_counts[stoi(split_line[3])] += 1;
            }
        }

        if (id == 0) {
            reverse(kmer_prg->nodes.begin(), kmer_prg->nodes.end());
        }

        id = 0;
        for (const auto &n : kmer_prg->nodes) {
            assert(kmer_prg->nodes[id]->id == id);
            id++;
            assert(n->id < outnode_counts.size() or assert_msg(n->id << ">=" << outnode_counts.size()));
            assert(n->id < innode_counts.size() or assert_msg(n->id << ">=" << innode_counts.size()));
            n->out_nodes.reserve(outnode_counts[n->id]);
            n->in_nodes.reserve(innode_counts[n->id]);
        }

        myfile.clear();
        myfile.seekg(0, myfile.beg);

        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                if (split_line[2] == split_line[4]) {
                    from = std::stoi(split_line[1]);
                    to = std::stoi(split_line[3]);
                } else {
                    // never happens
                    from = std::stoi(split_line[3]);
                    to = std::stoi(split_line[1]);
                }
                kmer_prg->add_edge(kmer_prg->nodes[from], kmer_prg->nodes[to]);
                //nodes[from]->out_nodes.push_back(nodes.at(to));
                //nodes[to]->in_nodes.push_back(nodes.at(from));
            }
        }
    } else {
        std::cerr << "Unable to open kmergraph file " << filepath << std::endl;
        exit(1);
    }
}