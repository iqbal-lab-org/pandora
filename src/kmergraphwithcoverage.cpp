#include <sstream>
#include <limits>
#include <cstdlib> /* srand, rand */
#include <cmath>

#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/log/trivial.hpp>

#include "kmergraphwithcoverage.h"
#include "localPRG.h"

using namespace prg;

void KmerGraphWithCoverage::set_exp_depth_covg(const uint32_t edp)
{
    const bool exp_depth_covg_parameter_is_valid = edp > 0;
    if (!exp_depth_covg_parameter_is_valid) {
        fatal_error(
            "Error setting exp_depth_covg: exp_depth_covg is invalid, must be > 0, is ",
            edp);
    }
    exp_depth_covg = edp;
}

void KmerGraphWithCoverage::set_binomial_parameter_p(const float e_rate)
{
    BOOST_LOG_TRIVIAL(debug) << "Set p in kmergraph";

    const bool valid_parameters_to_set_p
        = (kmer_prg->k != 0) && (0 < e_rate and e_rate < 1);
    if (!valid_parameters_to_set_p) {
        fatal_error("Error setting binomial parameter p, invalid parameters: ",
            "kmer_prg->k = ", kmer_prg->k, ", e_rate = ", e_rate);
    }

    binomial_parameter_p = 1 / exp(e_rate * kmer_prg->k);
}

void KmerGraphWithCoverage::increment_covg(
    uint32_t node_id, pandora::Strand strand, uint32_t sample_id)
{
    const bool sample_is_valid
        = this->node_index_to_sample_coverage[node_id].size() > sample_id;
    if (!sample_is_valid) {
        fatal_error(
            "Error incrementing coverage: sample_id is invalid (", sample_id, ")");
    }

    // get a pointer to the value we want to increment
    uint16_t* coverage_ptr = nullptr;
    if (strand == pandora::Strand::Forward) {
        coverage_ptr = &(this->node_index_to_sample_coverage[node_id][sample_id].first);
    } else {
        coverage_ptr
            = &(this->node_index_to_sample_coverage[node_id][sample_id].second);
    }

    const bool safe_to_increase_covg { (*coverage_ptr) < UINT16_MAX };
    if (safe_to_increase_covg) {
        ++(*coverage_ptr);
    }
}

uint32_t KmerGraphWithCoverage::get_covg(
    uint32_t node_id, pandora::Strand strand, uint32_t sample_id) const
{

    if (this->node_index_to_sample_coverage[node_id].size() <= sample_id)
        return 0;

    if (strand == pandora::Strand::Forward) {
        return (uint32_t)(
            this->node_index_to_sample_coverage[node_id][sample_id].first);
    } else {
        return (uint32_t)(
            this->node_index_to_sample_coverage[node_id][sample_id].second);
    }
}

void KmerGraphWithCoverage::set_covg(
    uint32_t node_id, uint16_t value, pandora::Strand strand, uint32_t sample_id)
{
    const bool sample_is_valid
        = this->node_index_to_sample_coverage[node_id].size() > sample_id;
    if (!sample_is_valid) {
        fatal_error("Error setting coverage: sample_id is invalid (", sample_id, ")");
    }

    if (strand == pandora::Strand::Forward) {
        this->node_index_to_sample_coverage[node_id][sample_id].first = value;
    } else {
        this->node_index_to_sample_coverage[node_id][sample_id].second = value;
    }
}

void KmerGraphWithCoverage::set_negative_binomial_parameters(
    const float& nbin_prob, const float& nb_fail)
{
    if (nbin_prob == 0 and nb_fail == 0)
        return;

    const bool negative_binomial_parameters_were_previously_set
        = (negative_binomial_parameter_p > 0 and negative_binomial_parameter_p < 1)
        && (negative_binomial_parameter_r > 0);
    if (!(negative_binomial_parameters_were_previously_set)) {
        fatal_error("Error setting negative_binomial_parameters: "
                    "negative_binomial_parameter_p (",
            negative_binomial_parameter_p, ")", " or negative_binomial_parameter_r (",
            negative_binomial_parameter_r, ") ", "were not correctly set");
    }

    negative_binomial_parameter_p += nbin_prob;
    negative_binomial_parameter_r += nb_fail;
}

float KmerGraphWithCoverage::nbin_prob(uint32_t node_id, const uint32_t& sample_id)
{
    const auto k = this->get_forward_covg(node_id, sample_id)
        + this->get_reverse_covg(node_id, sample_id);
    const float prob = pdf(boost::math::negative_binomial(
        negative_binomial_parameter_r, negative_binomial_parameter_p), k);
    const float log_prob = log(prob);
    const float return_prob = std::max(log_prob, std::numeric_limits<float>::lowest() / 1000);
    return return_prob;
}

float KmerGraphWithCoverage::lin_prob(uint32_t node_id, const uint32_t& sample_id)
{
    const bool reads_were_mapped_to_this_kmer_graph = num_reads != 0;
    if (!reads_were_mapped_to_this_kmer_graph) {
        fatal_error(
            "Impossible to compute lin_prob, no reads were mapped to this kmer graph");
    }
    auto k = this->get_forward_covg(node_id, sample_id)
        + this->get_reverse_covg(node_id, sample_id);
    return log(float(k) / num_reads);
}

float KmerGraphWithCoverage::bin_prob(uint32_t node_id, const uint32_t& sample_id)
{
    const bool reads_were_mapped_to_this_kmer_graph = num_reads != 0;
    if (!reads_were_mapped_to_this_kmer_graph) {
        fatal_error(
            "Impossible to compute bin_prob, no reads were mapped to this kmer graph");
    }
    return bin_prob(node_id, num_reads, sample_id);
}

float KmerGraphWithCoverage::bin_prob(
    const uint32_t& node_id, const uint32_t& num, const uint32_t& sample_id)
{
    const bool binomial_parameter_p_is_set_correctly = binomial_parameter_p != 1;
    if (!binomial_parameter_p_is_set_correctly) {
        fatal_error("Error when computing bin_prob: binomial_parameter_p (",
            binomial_parameter_p, ") is not correctly set");
    }

    const bool node_exists = node_id < kmer_prg->nodes.size();
    if (!node_exists) {
        fatal_error("Error when computing bin_prob: attempt to access inexistent node ",
            node_id);
    }

    uint32_t sum_coverages = this->get_forward_covg(node_id, sample_id)
        + this->get_reverse_covg(node_id, sample_id);

    float return_prob;
    if (node_id == (*(kmer_prg->sorted_nodes.begin()))->id
        or node_id == (*(kmer_prg->sorted_nodes.rbegin()))->id) {
        return_prob = 0; // is really undefined
    } else if (sum_coverages > num) {
        // under model assumptions this can't happen, but it inevitably will, so bodge
        return_prob
            = lognchoosek2(sum_coverages, this->get_forward_covg(node_id, sample_id),
                  this->get_reverse_covg(node_id, sample_id))
            + sum_coverages * log(binomial_parameter_p / 2);
        // note this may give disadvantage to repeat kmers
    } else {
        return_prob = lognchoosek2(num, this->get_forward_covg(node_id, sample_id),
                          this->get_reverse_covg(node_id, sample_id))
            + sum_coverages * log(binomial_parameter_p / 2)
            + (num - sum_coverages) * log(1 - binomial_parameter_p);
    }
    return return_prob;
}

float KmerGraphWithCoverage::get_prob(
    const std::string& prob_model, const uint32_t& node_id, const uint32_t& sample_id)
{
    if (prob_model == "nbin") {
        // is there no parameter check here?
        return nbin_prob(node_id, sample_id);
    } else if (prob_model == "bin") {
        const bool binomial_parameters_are_ok
            = (binomial_parameter_p < 1) && (num_reads > 0);
        if (!binomial_parameters_are_ok) {
            fatal_error("Error when computing kmer prob: binomial parameters are not "
                        "ok (binomial_parameter_p = ",
                binomial_parameter_p, ", ", "num_reads = ", num_reads);
        }
        return bin_prob(node_id, sample_id);
    } else if (prob_model == "lin") {
        // is there no parameter check here?
        return lin_prob(node_id, sample_id);
    } else {
        fatal_error("Invalid probability model for kmer coverage distribution: ",
            prob_model, ". Should be nbin, bin or lin");
    }
}

bool KmerGraphWithCoverage::coverage_is_zeroes(const uint32_t& sample_id)
{
    bool all_zero = true;
    for (const auto& node_ptr : kmer_prg->nodes) {
        const auto covg { this->get_forward_covg(node_ptr->id, sample_id)
            + this->get_reverse_covg(node_ptr->id, sample_id) };
        if (covg > 0) {
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

float KmerGraphWithCoverage::find_max_path(std::vector<KmerNodePtr>& maxpath,
    const std::string& prob_model, const uint32_t& max_num_kmers_to_average,
    const uint32_t& sample_id, const pangenome::Node *pangenome_node)
{
    // TODO: FIX THIS INNEFICIENCY I INTRODUCED
    const std::vector<KmerNodePtr> sorted_nodes(
        this->kmer_prg->sorted_nodes.begin(), this->kmer_prg->sorted_nodes.end());
    this->kmer_prg->check();

    // also check not all 0 covgs
    auto coverages_all_zero = coverage_is_zeroes(sample_id);
    if (coverages_all_zero)
        return std::numeric_limits<float>::lowest();

    // create vectors to hold the intermediate values
    std::vector<float> max_sum_of_log_probs_from_node(sorted_nodes.size(), 0);
    std::vector<uint32_t> length_of_maxpath_from_node(sorted_nodes.size(), 0);
    std::vector<uint32_t> prev_node_along_maxpath(
        sorted_nodes.size(), sorted_nodes.size() - 1);
    float max_mean;
    float max_sum;
    int max_length;
    const float float_point_tolerance = 0.000001;
    const float likelihood_boost_per_node = 0.005;  // if a path is longer, we give a 0.5% log likelihood boost for each additional node
    const float likelihood_boost_ceiling = 0.50; // the maximum log likelihood boost is 50%

    for (uint32_t j = sorted_nodes.size() - 1; j != 0; --j) {
        max_mean = std::numeric_limits<float>::lowest();
        max_sum = std::numeric_limits<float>::lowest();
        max_length = 1; // tie break with longest kmer path
        const auto& current_node = sorted_nodes[j - 1];
        if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " current_node " << pangenome_node->prg->string_along_path(current_node->path);
        for (uint32_t i = 0; i != current_node->out_nodes.size(); ++i) {
            const auto& considered_outnode = current_node->out_nodes[i].lock();
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode " << pangenome_node->prg->string_along_path(considered_outnode->path);

            const bool is_terminus_and_most_likely
                = considered_outnode->id == sorted_nodes.back()->id
                and thresh > max_mean + float_point_tolerance;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is_terminus_and_most_likely? " << is_terminus_and_most_likely;

            float considered_outnode_mean = std::numeric_limits<float>::lowest();
            if (length_of_maxpath_from_node[considered_outnode->id])
                considered_outnode_mean = max_sum_of_log_probs_from_node[considered_outnode->id] / length_of_maxpath_from_node[considered_outnode->id];

            const int delta_length = int(length_of_maxpath_from_node[considered_outnode->id]) - max_length;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode delta_length " << delta_length;

            const bool is_longer_path = delta_length > 0;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is_longer_path? " << is_longer_path;
            const bool is_shorter_path = delta_length < 0;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is_shorter_path? " << is_shorter_path;

            float considered_outnode_mean_with_boost = considered_outnode_mean;
            float max_mean_with_boost = max_mean;
            float likelihood_boost = delta_length * likelihood_boost_per_node;
            likelihood_boost = std::min(likelihood_boost, likelihood_boost_ceiling);
            if (is_longer_path) {
                considered_outnode_mean_with_boost = considered_outnode_mean *
                    (1 - likelihood_boost);
            }
            if (is_shorter_path) {
                max_mean_with_boost = max_mean * (1 - likelihood_boost);
            }

            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode considered_outnode_mean " << considered_outnode_mean;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode considered_outnode_mean_with_boost " << considered_outnode_mean_with_boost;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode max_mean " << max_mean;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode max_mean_with_boost " << max_mean_with_boost;

            const bool avg_log_likelihood_is_most_likely = considered_outnode_mean_with_boost > max_mean_with_boost + float_point_tolerance;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode avg_log_likelihood_is_most_likely? " << avg_log_likelihood_is_most_likely;

            const bool avg_log_likelihood_is_most_likely_without_boost =
                considered_outnode_mean > max_mean + float_point_tolerance;
            if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode avg_log_likelihood_is_most_likely_without_boost? " << avg_log_likelihood_is_most_likely_without_boost;

            if (avg_log_likelihood_is_most_likely != avg_log_likelihood_is_most_likely_without_boost) {
                if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode avg_log_likelihood_is_most_likely != avg_log_likelihood_is_most_likely_without_boost";
            }

            if (is_terminus_and_most_likely or avg_log_likelihood_is_most_likely) {
                if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is_terminus_and_most_likely or avg_log_likelihood_is_most_likely";
                float current_node_prob = get_prob(prob_model, current_node->id, sample_id);
                current_node_prob = current_node_prob * likelihood_boost_per_node;
                max_sum_of_log_probs_from_node[current_node->id] = current_node_prob
                    + max_sum_of_log_probs_from_node[considered_outnode->id];
                length_of_maxpath_from_node[current_node->id]
                    = 1 + length_of_maxpath_from_node[considered_outnode->id];
                prev_node_along_maxpath[current_node->id] = considered_outnode->id;

                /*
                if (length_of_maxpath_from_node[current_node->id]
                    > max_num_kmers_to_average) {
                    uint32_t prev_node = prev_node_along_maxpath[current_node->id];
                    for (uint step = 0; step < max_num_kmers_to_average; step++) {
                        prev_node = prev_node_along_maxpath[prev_node];
                    }
                    max_sum_of_log_probs_from_node[current_node->id]
                        -= get_prob(prob_model, sorted_nodes[prev_node]->id, sample_id);
                    length_of_maxpath_from_node[current_node->id] -= 1;

                    // this remains as an assert, as it is a code check
                    // Note: I think we might even be able to remove this
                    assert(length_of_maxpath_from_node[current_node->id]
                        == max_num_kmers_to_average);
                }
                 */

                if (considered_outnode->id != sorted_nodes.back()->id) {
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is not terminus";
                    max_mean = considered_outnode_mean;
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode new max_mean " << max_mean;
                    max_sum = max_sum_of_log_probs_from_node[considered_outnode->id];
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode new max_sum " << max_sum;
                    max_length = length_of_maxpath_from_node[considered_outnode->id];
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode new max_length " << max_length;
                } else {
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode is terminus";
                    max_mean = thresh;
                    if (pangenome_node) BOOST_LOG_TRIVIAL(debug) << "[Sample " << sample_id << " ML_path_algorithm] " << pangenome_node->name << " considered_outnode new max_mean " << max_mean;
                }
            }
        }
    }

    // extract path
    uint32_t prev_node = prev_node_along_maxpath[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1) {
        maxpath.push_back(this->kmer_prg->nodes[prev_node]);
        prev_node = prev_node_along_maxpath[prev_node];

        if (maxpath.size() > 1000000) {
            fatal_error("I think I've found an infinite loop - is "
                        "something wrong with this kmergraph?");
        }
    }

    const bool path_was_found_through_the_kmer_PRG = length_of_maxpath_from_node[0] > 0;
    if (!path_was_found_through_the_kmer_PRG) {
        fatal_error("Error when finding max path: found no path through kmer prg");
    }

    return prob_path(maxpath, sample_id, prob_model);
}

std::vector<std::vector<KmerNodePtr>> KmerGraphWithCoverage::get_random_paths(
    uint32_t num_paths)
{
    // find a random path through kmergraph picking ~uniformly from the outnodes at each
    // point
    std::vector<std::vector<KmerNodePtr>> rpaths;
    std::vector<KmerNodePtr> rpath;
    uint32_t i;

    time_t now;
    now = time(nullptr);
    srand((unsigned int)now);

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

float KmerGraphWithCoverage::prob_path(const std::vector<KmerNodePtr>& kpath,
    const uint32_t& sample_id, const std::string& prob_model)
{
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

void KmerGraphWithCoverage::save_covg_dist(const std::string& filepath)
{
    std::ofstream handle;
    handle.open(filepath);

    for (const auto& kmer_node_ptr : kmer_prg->nodes) {
        const KmerNode& kmer_node = *kmer_node_ptr;

        uint32_t sample_id = 0;
        for (const auto& sample_coverage :
            node_index_to_sample_coverage[kmer_node.id]) {
            handle << kmer_node.id << " " << sample_id << " " << sample_coverage.first
                   << " " << sample_coverage.second;

            sample_id++;
        }
    }
    handle.close();
}

// save the KmerGraph as gfa
// TODO: THIS SHOULD BE RECODED, WE ARE DUPLICATING CODE HERE (SEE KmerGraph::save())!!!
void KmerGraphWithCoverage::save(
    const fs::path& filepath, const std::shared_ptr<LocalPRG> localprg) const
{
    uint32_t sample_id = 0;

    fs::ofstream handle(filepath);
    if (handle.is_open()) {
        handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
        for (const auto& c : kmer_prg->nodes) {
            handle << "S\t" << c->id << "\t";

            if (localprg != nullptr) {
                handle << localprg->string_along_path(c->path);
            } else {
                handle << c->path;
            }

            handle << "\tFC:i:" << this->get_forward_covg(c->id, sample_id)
                   << "\tRC:i:" << this->get_reverse_covg(c->id, sample_id)
                   << std::endl;

            for (uint32_t j = 0; j < c->out_nodes.size(); ++j) {
                handle << "L\t" << c->id << "\t+\t" << c->out_nodes[j].lock()->id
                       << "\t+\t0M" << std::endl;
            }
        }
        handle.close();
    } else {
        fatal_error("Unable to open kmergraph file ", filepath);
    }
}

// TODO: THIS SHOULD BE RECODED, WE ARE DUPLICATING CODE HERE (SEE KmerGraph::load())!!!
// TODO: remove this method?
void KmerGraphWithCoverage::load(const std::string& filepath)
{
    // TODO: this might be dangerous, recode this?
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

                const bool line_is_consistent = split_line.size() >= 4;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                id = std::stoi(split_line[1]);
                num_nodes = std::max(num_nodes, id);
            }
        }
        myfile.clear();
        myfile.seekg(0, myfile.beg);
        kmer_prg->nodes.reserve(num_nodes);
        std::vector<uint16_t> outnode_counts(num_nodes + 1, 0),
            innode_counts(num_nodes + 1, 0);

        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");

                const bool line_is_consistent = split_line.size() >= 4;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                id = stoi(split_line[1]);
                ss << split_line[2];
                char c = ss.peek();

                if (!isdigit(c)) {
                    fatal_error("Error reading GFA: cannot read in this sort of "
                                "kmergraph GFA as it ",
                        "does not label nodes with their PRG path. ",
                        "Offending line: ", line);
                }

                ss >> p;
                ss.clear();
                // add_node(p);
                KmerNodePtr n = std::make_shared<KmerNode>(id, p);

                const bool id_is_consistent = (id == kmer_prg->nodes.size()
                    or num_nodes - id == kmer_prg->nodes.size());
                if (!id_is_consistent) {
                    fatal_error("Error reading GFA: node ID is inconsistent.",
                        "id = ", id, ", ", "nodes.size() = ", kmer_prg->nodes.size(),
                        ", ", "num_nodes = ", num_nodes);
                }

                kmer_prg->nodes.push_back(n);
                kmer_prg->sorted_nodes.insert(n);
                if (kmer_prg->k == 0 and p.length() > 0) {
                    kmer_prg->k = p.length();
                }
                covg = (uint16_t)(stoul(split(split_line[3], "FC:i:")[0]));
                this->set_forward_covg(n->id, covg, sample_id);
                covg = (uint16_t)(stoul(split(split_line[4], "RC:i:")[0]));
                this->set_reverse_covg(n->id, covg, sample_id);
                if (split_line.size() >= 6) {
                    n->num_AT = std::stoi(split_line[5]);
                }
            } else if (line[0] == 'L') {
                split_line = split(line, "\t");

                const bool line_is_consistent = split_line.size() >= 5;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                const int from_node = stoi(split_line[1]);
                const int to_node = stoi(split_line[3]);
                const bool from_node_in_range = from_node < (int)outnode_counts.size();
                const bool to_node_in_range = to_node < (int)innode_counts.size();
                if (!from_node_in_range) {
                    fatal_error(
                        "Error reading GFA: from_node out of range: ", from_node,
                        ">=", outnode_counts.size(), ". Offending line: ", line);
                }
                if (!to_node_in_range) {
                    fatal_error("Error reading GFA: to_node out of range: ", to_node,
                        ">=", innode_counts.size(), ". Offending line: ", line);
                }

                outnode_counts[stoi(split_line[1])] += 1;
                innode_counts[stoi(split_line[3])] += 1;
            }
        }

        if (id == 0) {
            reverse(kmer_prg->nodes.begin(), kmer_prg->nodes.end());
        }

        id = 0;
        for (const auto& n : kmer_prg->nodes) {
            const bool id_is_consistent = (kmer_prg->nodes[id]->id == id)
                && (n->id < outnode_counts.size()) && (n->id < innode_counts.size());
            if (!id_is_consistent) {
                fatal_error("Error reading GFA: node: ", n,
                    " has inconsistent id, should be ", id);
            }
            id++;
            n->out_nodes.reserve(outnode_counts[n->id]);
            n->in_nodes.reserve(innode_counts[n->id]);
        }

        myfile.clear();
        myfile.seekg(0, myfile.beg);

        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");

                const bool line_is_consistent = split_line.size() >= 5;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                if (split_line[2] == split_line[4]) {
                    from = std::stoi(split_line[1]);
                    to = std::stoi(split_line[3]);
                } else {
                    // never happens
                    from = std::stoi(split_line[3]);
                    to = std::stoi(split_line[1]);
                }
                kmer_prg->add_edge(kmer_prg->nodes[from], kmer_prg->nodes[to]);
            }
        }
    } else {
        fatal_error("Error reading GFA: unable to open kmergraph file: ", filepath);
    }
}