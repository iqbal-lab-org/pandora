#include <iomanip>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <cmath>
#include <set>
#include <memory>
#include <ctime>
#include <algorithm>
#include <boost/filesystem.hpp>

#include "utils.h"
#include "seq.h"
#include "localPRG.h"
#include "pangenome/pangraph.h"
#include "minihit.h"
#include "fastaq_handler.h"
#include "minimizermatch_file.h"
#include "sam_file.h"
#include "minihit_clusters.h"

#define _GNU_SOURCE
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <array>

std::string now()
{
    time_t now;
    std::string dt;

    now = time(nullptr);
    dt = ctime(&now);
    return dt.substr(0, dt.length() - 1) + " ";
}

std::string int_to_string(const int number)
{
    std::stringstream ss;
    ss << std::setw(2) << std::setfill('0') << number;
    return ss.str();
}

std::vector<std::string> split(const std::string& query, const std::string& d)
{
    std::vector<std::string> v;
    std::string::size_type k = 0;
    std::string::size_type j = query.find(d, k);
    while (j != std::string::npos) {
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

char complement(char n)
{
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
    return 'N';
}

std::string rev_complement(std::string s)
{
    transform(begin(s), end(s), begin(s), complement);
    reverse(s.begin(), s.end());
    return s;
}

float lognchoosek2(uint32_t n, uint32_t k1, uint32_t k2)
{
    const bool parameters_are_valid = n >= (k1 + k2);
    if (!parameters_are_valid) {
        fatal_error(
            "Currently the model assumes that the most a given kmer (defined by "
            "position) can occur is once per read, i.e. an error somewhere else in the "
            "read cannot result in this kmer. If you are getting this message, then "
            "you have evidence of violation of this assumption. Either try using a "
            "bigger k, or come up with a better model");
    }
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

void load_vcf_refs_file(const fs::path& filepath, VCFRefs& vcf_refs)
{
    BOOST_LOG_TRIVIAL(info) << "Loading VCF refs from file " << filepath;

    FastaqHandler fh(filepath.string());
    while (!fh.eof()) {
        try {
            fh.get_next();
        } catch (std::out_of_range& err) {
            break;
        }
        if (!fh.name.empty() && !fh.read.empty()) {
            vcf_refs[fh.name] = fh.read;
        }
    }
}

void add_read_hits(const Seq& sequence,
    const std::shared_ptr<MinimizerHits>& minimizer_hits, const Index& index)
{
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and
    // adds hits to a global MinimizerHits object
    // Seq s(id, name, seq, w, k);
    for (auto sequenceSketchIt = sequence.sketch.begin();
         sequenceSketchIt != sequence.sketch.end(); ++sequenceSketchIt) {
        auto minhashIt = index.minhash.find((*sequenceSketchIt).canonical_kmer_hash);
        if (minhashIt != index.minhash.end()) { // checks if the kmer is in the index
            // yes, add all hits of this minimizer hit to this kmer
            for (const MiniRecord& miniRecord : *(minhashIt->second)) {
                minimizer_hits->insert(sequence.id, *sequenceSketchIt, miniRecord);
            }
        }
    }
}


void decide_if_add_cluster_or_not(
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits, // Note: clusters_of_hits here is in insertion mode
    std::vector<uint32_t> &prg_max_path_lengths,
    const std::vector<std::string> &prg_names,
    const std::set<MinimizerHitPtr, pComp>::iterator &mh_previous,
    const uint32_t expected_number_kmers_in_read_sketch,
    const float fraction_kmers_required_for_cluster,
    const uint32_t min_cluster_size,
    const MinimizerHits &current_cluster,
    const std::vector<uint32_t> &distances_between_hits,
    ClusterDefFile& cluster_def_file) {
    // keep clusters which cover at least 1/2 the expected number of minihits
    const uint32_t length_based_threshold =
        std::min(prg_max_path_lengths[(*mh_previous)->get_prg_id()], expected_number_kmers_in_read_sketch)
                                      * fraction_kmers_required_for_cluster;

    const uint32_t cluster_size_threshold = std::max(length_based_threshold, min_cluster_size);
    const bool cluster_should_be_accepted = current_cluster.get_number_of_unique_mini_in_cluster() >= cluster_size_threshold;

    if (cluster_should_be_accepted) {
        clusters_of_hits.insert(current_cluster);
    }

    // Note: this is a slow critical region and could be optimised, but there is no need
    // to, as this is just run when debugging files should be created, and is expected
    // to be slow.
    if (!cluster_def_file.is_fake_file) {
#pragma omp critical(cluster_def_file)
        {
            cluster_def_file << seq.name << "\t" << prg_names[(*mh_previous)->get_prg_id()] << "\t";
            if (cluster_should_be_accepted) {
                cluster_def_file << "accepted\t";
            }else {
                cluster_def_file << "rejected\t";
            }
            cluster_def_file << current_cluster.size() << "\t"
                             << current_cluster.get_number_of_equal_read_minimizers() << "\t"
                             << current_cluster.get_number_of_unique_mini_in_cluster() << "\t"
                             << length_based_threshold << "\t"
                             << min_cluster_size << "\t";

            for (auto current_cluster_it=current_cluster.begin();
                 current_cluster_it != current_cluster.end();
                 current_cluster_it++) {
                cluster_def_file << (*current_cluster_it)->get_read_start_position() << ",";
            }
            cluster_def_file << "\t";

            for (const auto &distance : distances_between_hits) {
                cluster_def_file << distance << ",";
            }
            cluster_def_file << "\n";
        }
    }
}

void define_clusters(
    const std::string &sample_name,
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits, // Note: clusters_of_hits here is in insertion mode
    std::vector<uint32_t> &prg_max_path_lengths,
    const std::vector<std::string> &prg_names,
    std::shared_ptr<MinimizerHits> &minimizer_hits, const int max_diff,
    const float& fraction_kmers_required_for_cluster, const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch,
    ClusterDefFile& cluster_def_file)
{
    const std::string tag = "[Sample: " + sample_name + ", read index: " + to_string(seq.id) + "]: ";

    BOOST_LOG_TRIVIAL(trace) << tag << "Define clusters of hits from the "
                             << minimizer_hits->size() << " hits";

    if (minimizer_hits->empty()) {
        return;
    }

    // A cluster of hits should match same localPRG, each hit not more than max_diff
    // read bases from the last hit (this last bit is to handle repeat genes).
    auto mh_previous = minimizer_hits->begin();
    MinimizerHits current_cluster(&prg_max_path_lengths);
    current_cluster.insert(*mh_previous);

    std::vector<uint32_t> distances_between_hits;
    for (auto mh_current = ++minimizer_hits->begin();
         mh_current != minimizer_hits->end(); ++mh_current) {
        const bool read_minimizer_is_the_same =
            (*mh_current)->get_read_start_position() == (*mh_previous)->get_read_start_position();
        if (read_minimizer_is_the_same)
            current_cluster.increment_number_of_equal_read_minimizers();

        const bool switched_reads = (*mh_current)->get_read_id() != (*mh_previous)->get_read_id();
        const bool switched_prgs = (*mh_current)->get_prg_id() != (*mh_previous)->get_prg_id();
        const uint32_t distance_between_hits = abs((int)(*mh_current)->get_read_start_position()
                                                   - (int)(*mh_previous)->get_read_start_position());
        const bool hits_too_distant = distance_between_hits > max_diff;

        const bool unconsistent_strands
            = (*mh_current)->same_strands() != (*mh_previous)->same_strands();

        const bool switched_clusters = switched_reads or switched_prgs or
            hits_too_distant or unconsistent_strands;
        if (switched_clusters) {
            decide_if_add_cluster_or_not(seq, clusters_of_hits, prg_max_path_lengths, prg_names,
                mh_previous, expected_number_kmers_in_read_sketch,
                fraction_kmers_required_for_cluster, min_cluster_size, current_cluster,
                distances_between_hits,cluster_def_file);

            // prepare next cluster
            current_cluster.clear();
            distances_between_hits.clear();

        } else {
            if (!read_minimizer_is_the_same) {
                distances_between_hits.push_back(distance_between_hits);
            }
        }
        current_cluster.insert(*mh_current);
        mh_previous = mh_current;
    }
    decide_if_add_cluster_or_not(seq, clusters_of_hits, prg_max_path_lengths, prg_names, mh_previous,
                                 expected_number_kmers_in_read_sketch, fraction_kmers_required_for_cluster,
                                 min_cluster_size, current_cluster,
                                 distances_between_hits, cluster_def_file);
}

MinimizerHitClusters filter_clusters(
    const std::string &sample_name,
    const Seq &seq,
    const MinimizerHitClusters& clusters_of_hits,
    const std::vector<std::string> &prg_names,
    ClusterFilterFile& cluster_filter_file,
    const double overlap_threshold,
    const uint32_t rng_seed)
{
    const std::string tag = "[Sample: " + sample_name + ", read index: " + to_string(seq.id) + "]: ";
    MinimizerHitClusters filtered_clusters_of_hits(rng_seed);

    // Next order clusters, choose between those that overlap by too much
    BOOST_LOG_TRIVIAL(trace) << tag << "Filter the " << clusters_of_hits.size()
                             << " clusters of hits";
    if (clusters_of_hits.empty()) {
        filtered_clusters_of_hits.finalise_insertions();
        return filtered_clusters_of_hits;
    }

    std::set<size_t> clusters_to_remove;

    // to do this consider pairs of clusters in turn
    auto c_previous = clusters_of_hits.begin();
    for (auto c_current = ++clusters_of_hits.begin();
         c_current != clusters_of_hits.end(); ++c_current) {
        const bool both_clusters_from_same_read = c_current->front()->get_read_id() == c_previous->front()->get_read_id();
        const bool same_prg = c_current->front()->get_prg_id() == c_previous->front()->get_prg_id();
        const bool unconsistent_strands = c_current->front()->same_strands() != c_previous->front()->same_strands();
        const double overlap = c_current->overlap_amount(*c_previous);
        const bool clusters_overlap = overlap >= overlap_threshold;
        const bool one_of_the_clusters_should_be_filtered_out =
            both_clusters_from_same_read && ((same_prg && unconsistent_strands) or
                clusters_overlap);
        if (one_of_the_clusters_should_be_filtered_out)
        // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters,
        // but could also impose no more than 2k hits in overlap
        {
            const bool should_remove_current_cluster = c_previous->is_preferred_to(*c_current);

            // Note: this is a slow critical region and could be optimised, but there is no need
            // to, as this is just run when debugging files should be created, and is expected
            // to be slow.
            if (!cluster_filter_file.is_fake_file) {
#pragma omp critical(cluster_filter_file)
                {
                    if (should_remove_current_cluster) {
                        cluster_filter_file
                            << seq.name << "\t"
                            << prg_names[c_current->front()->get_prg_id()] << "\t"
                            << c_current->get_number_of_unique_mini_in_cluster() << "\t"
                            << c_current->target_coverage() << "\t"
                            << c_current->front()->get_read_start_position() << "\t"
                            << c_current->back()->get_read_start_position() << "\t"
                            << overlap << "\t"
                            << "filtered_out\t"
                            << prg_names[c_previous->front()->get_prg_id()] << "\t"
                            << c_previous->get_number_of_unique_mini_in_cluster() << "\t"
                            << c_previous->target_coverage() << "\t"
                            << c_previous->front()->get_read_start_position() << "\t"
                            << c_previous->back()->get_read_start_position() << "\n";
                    } else {
                        cluster_filter_file
                            << seq.name << "\t"
                            << prg_names[c_previous->front()->get_prg_id()] << "\t"
                            << c_previous->get_number_of_unique_mini_in_cluster() << "\t"
                            << c_previous->target_coverage() << "\t"
                            << c_previous->front()->get_read_start_position() << "\t"
                            << c_previous->back()->get_read_start_position() << "\t"
                            << overlap << "\t"
                            << "kept\t"
                            << prg_names[c_current->front()->get_prg_id()] << "\t"
                            << c_current->get_number_of_unique_mini_in_cluster() << "\t"
                            << c_current->target_coverage() << "\t"
                            << c_current->front()->get_read_start_position() << "\t"
                            << c_current->back()->get_read_start_position() << "\n";
                    }
                }
            }

            if (should_remove_current_cluster) {
                auto pos = c_current - clusters_of_hits.begin();
                clusters_to_remove.insert(pos);
                BOOST_LOG_TRIVIAL(trace) << tag << "Cluster #" << pos << " to be filtered out";
                // c_previous continues the same
            } else {
                auto pos = c_previous - clusters_of_hits.begin();
                clusters_to_remove.insert(pos);
                c_previous = c_current;
                BOOST_LOG_TRIVIAL(trace) << tag << "Cluster #" << pos << " to be filtered out";
            }
        }
        else {
            c_previous = c_current;
        }
    }

    size_t cluster_index = 0;
    for (const auto &cluster : clusters_of_hits) {
        const bool cluster_should_be_added = clusters_to_remove.find(cluster_index) == clusters_to_remove.end();
        if (cluster_should_be_added) {
            filtered_clusters_of_hits.insert(cluster);
        }
        ++cluster_index;
    }
    filtered_clusters_of_hits.finalise_insertions();

    // Note: this is a slow critical region and could be optimised, but there is no need
    // to, as this is just run when debugging files should be created, and is expected
    // to be slow.
    if (!cluster_filter_file.is_fake_file) {
#pragma omp critical(cluster_filter_file)
        {
            for (const auto& cluster : filtered_clusters_of_hits) {
                cluster_filter_file << seq.name << "\t"
                                    << prg_names[cluster.front()->get_prg_id()] << "\t"
                                    << cluster.get_number_of_unique_mini_in_cluster() << "\t"
                                    << cluster.target_coverage() << "\t"
                                    << cluster.front()->get_read_start_position() << "\t"
                                    << cluster.back()->get_read_start_position() << "\t"
                                    << "NA\t"
                                    << "kept\n";
            }
        }
    }

    BOOST_LOG_TRIVIAL(trace) << tag << "Now have " << filtered_clusters_of_hits.size()
                             << " clusters of hits";

    return filtered_clusters_of_hits;
}

void filter_clusters2(MinimizerHitClusters& clusters_of_hits,
    const uint32_t& genome_size)
{
    // TODO: this method is all commented out, tagging it for removal
    // Currently let's just error out if we ever call it
    fatal_error("Not implemented");
    /*
    // Sort clusters by size, and filter out those small clusters which are entirely
    // contained in bigger clusters on reads
    BOOST_LOG_TRIVIAL(trace) << "Filter2 the " << clusters_of_hits.size()
                             << " clusters of hits";
    if (clusters_of_hits.empty()) {
        return;
    }

    MinimizerHitClusters clusters_by_size(
        clusters_of_hits.begin(), clusters_of_hits.end());

    auto it = clusters_by_size.begin();
    std::vector<int> read_v(genome_size, 0);
    fill(read_v.begin() + (*(it->begin()))->get_read_start_position(),
        read_v.begin() + (*--(it->end()))->get_read_start_position(), 1);
    bool contained;
    for (auto it_next = ++clusters_by_size.begin(); it_next != clusters_by_size.end();
         ++it_next) {
        if ((*(it_next->begin()))->get_read_id() == (*(it->begin()))->get_read_id()) {
            // check if have any 0s in interval of read_v between first and last
            contained = true;
            for (uint32_t i = (*(it_next->begin()))->get_read_start_position();
                 i < (*--(it_next->end()))->get_read_start_position(); ++i) {
                if (read_v[i] == 0) {
                    contained = false;
                    fill(read_v.begin() + i,
                        read_v.begin()
                            + (*--(it_next->end()))->get_read_start_position(),
                        1);
                    break;
                }
            }
            if (contained) {
                clusters_of_hits.erase(*it_next);
            }
        } else {
            // consider new read
            fill(read_v.begin(), read_v.end(), 0);
        }
        ++it;
    }
    BOOST_LOG_TRIVIAL(trace) << "Now have " << clusters_of_hits.size()
                             << " clusters of hits";
    */
}

void add_clusters_to_pangraph(
    const MinimizerHitClusters& minimizer_hit_clusters,
    std::shared_ptr<pangenome::Graph> &pangraph,
    Index &index, uint32_t sample_id)
{
    BOOST_LOG_TRIVIAL(trace) << "Add clusters to PanGraph";
    if (minimizer_hit_clusters.empty()) {
        return;
    }

    for (const auto &cluster : minimizer_hit_clusters) {
        // each cluster here defines a mapping, so we know which prgs mapped
        // we lazily load them just now
        uint32_t mapped_prg_id = (*cluster.begin())->get_prg_id();
        std::shared_ptr<LocalPRG> prg = index.get_prg_given_id(mapped_prg_id);
        pangraph->record_hit(prg);
        auto& pangraph_node = pangraph->get_node(prg);
        for (const auto &hit : cluster) {
            if (hit->same_strands()) {
                pangraph_node->kmer_prg_with_coverage.increment_forward_covg(
                    hit->get_kmer_node_id(), sample_id);
            } else {
                pangraph_node->kmer_prg_with_coverage.increment_reverse_covg(
                    hit->get_kmer_node_id(), sample_id);
            }
        }
    }
}

MinimizerHitClusters get_minimizer_hit_clusters(
    const std::string &sample_name,
    const Seq &seq,
    std::vector<uint32_t> &prg_max_path_lengths,
    const std::vector<std::string> &prg_names,
    std::shared_ptr<MinimizerHits> &minimizer_hits,
    const int max_diff,
    const float& fraction_kmers_required_for_cluster,
    ClusterDefFile &cluster_def_file,
    ClusterFilterFile &cluster_filter_file,
    const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch,
    const uint32_t rng_seed)
{
    const std::string tag = "[Sample: " + sample_name + ", read index: " + to_string(seq.id) + "]: ";
    MinimizerHitClusters minimizer_hit_clusters(rng_seed);

    if (minimizer_hits->empty()) {
        minimizer_hit_clusters.finalise_insertions();
        BOOST_LOG_TRIVIAL(trace) << tag << "Found 0 clusters of hits";
        return minimizer_hit_clusters;
    }

    define_clusters(sample_name, seq, minimizer_hit_clusters, prg_max_path_lengths,
        prg_names, minimizer_hits, max_diff, fraction_kmers_required_for_cluster,
        min_cluster_size, expected_number_kmers_in_read_sketch, cluster_def_file);

    minimizer_hit_clusters.finalise_insertions();
    BOOST_LOG_TRIVIAL(trace) << tag << "Found " << minimizer_hit_clusters.size() << " clusters of hits";

    MinimizerHitClusters filtered_clusters_of_hits = filter_clusters(sample_name, seq, minimizer_hit_clusters, prg_names, cluster_filter_file);
    // filter_clusters2(clusters_of_hits, genome_size);

    return filtered_clusters_of_hits;
}

// TODO: this should be in a constructor of pangenome::Graph or in a factory class
uint32_t pangraph_from_read_file(const SampleData& sample,
    std::shared_ptr<pangenome::Graph> &pangraph, Index &index,
    const int max_diff, const float& e_rate,
    const fs::path& sample_outdir, const uint32_t min_cluster_size,
    const uint32_t genome_size, const uint32_t max_covg, uint32_t threads,
    const bool keep_extra_debugging_files, const uint32_t rng_seed)
{
    // constant variables
    const SampleIdText sample_name = sample.first;
    const SampleFpath sample_filepath = sample.second;
    const std::string tag = "[Sample " + sample_name + "]: ";
    const uint32_t w = index.get_window_size();
    const uint32_t k = index.get_kmer_size();
    const double fraction_kmers_required_for_cluster = 0.5 / exp(e_rate * k);
    const uint32_t nb_reads_to_map_in_a_batch = 1000;

    BOOST_LOG_TRIVIAL(trace) << tag << "e_rate: " << e_rate;
    BOOST_LOG_TRIVIAL(trace) << tag << "k: " << k;
    BOOST_LOG_TRIVIAL(trace) << tag << "exp(e_rate * k): " << exp(e_rate * k);
    BOOST_LOG_TRIVIAL(trace) << tag << "fraction_kmers_required_for_cluster: " << fraction_kmers_required_for_cluster;

    // shared variable - controlled by critical(covg)
    uint64_t covg { 0 };

    // shared variables - controlled by critical(ReadFileMutex)
    FastaqHandler fh(sample_filepath);
    uint32_t id { 0 };

    SAMFile filtered_mappings(sample_outdir / (sample_name + ".filtered.sam"),
        index.get_prg_names(), index.get_prg_lengths(), k*2, k);

    MinimizerMatchFile minimizer_matches(sample_outdir / (sample_name + ".minimatches"),
        index.get_prg_names(), !keep_extra_debugging_files);
    PafFile paf_file(sample_outdir / (sample_name + ".minipaf"),
        index.get_prg_names(), !keep_extra_debugging_files);

    ClusterDefFile cluster_def_file(sample_outdir / (sample_name + ".clusters_def_report"), !keep_extra_debugging_files);
    ClusterFilterFile cluster_filter_file(sample_outdir / (sample_name + ".clusters_filter_report"), !keep_extra_debugging_files);


// parallel region
#pragma omp parallel num_threads(threads)
    {
        // will hold the reads batch
        std::vector<Seq> sequencesBuffer(
            nb_reads_to_map_in_a_batch, Seq(0, "null", "", w, k));
        while (true) {
            // read the next batch of reads
            uint32_t nbOfReads = 0;

// read the reads in batch
#pragma omp critical(ReadFileMutex)
            {
                for (auto& sequence : sequencesBuffer) {
                    if (id && id % 100000 == 0) {
                        BOOST_LOG_TRIVIAL(info) << id << " reads processed...";
                    }
                    try {
                        fh.get_next();
                    } catch (std::out_of_range& err) {
                        break;
                    }
                    sequence.initialize(id, fh.name, fh.read, w, k);
                    ++nbOfReads;
                    ++id;
                }
            }

            if (nbOfReads == 0)
                break; // we reached the end of the file, nothing else to map

            // quasimap the batch of reads
            bool coverageExceeded = false;
            for (uint32_t i = 0; i < nbOfReads; i++) {
                const auto& sequence = sequencesBuffer[i];

                // checks if we are still good regarding coverage
                if (!sequence.sketch.empty()) {
#pragma omp critical(covg)
                    {
                        // check if the max_covg was already exceeded
                        if (covg / genome_size > max_covg) {
                            // if reached here, it means that another thread realised
                            // that we went past the max_covg, so we just exit
                            coverageExceeded = true;
                        } else {
                            // no other thread still signalized exceeding max coverage
                            covg += sequence.length();
                            if (covg / genome_size > max_covg) {
                                // oops, we are the first one to see max_covg being
                                // exceeded, print and exit!
                                BOOST_LOG_TRIVIAL(warning)
                                    << "Stop processing reads as have reached max "
                                       "coverage";
                                coverageExceeded = true;
                            }
                        }
                    }

                    if (coverageExceeded)
                        break; // max covg exceeded, get out
                } else {
                    continue;
                }

                const auto expected_number_kmers_in_read_sketch { sequence.length() * 2
                    / (w + 1) };

                // get the minimizer hits
                auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits(&index.get_prg_max_path_lengths()));
                add_read_hits(sequence, minimizer_hits, index);

                // write unfiltered minimizer hits
                if (!minimizer_matches.is_fake_file) {
#pragma omp critical(minimizer_matches)
                    {
                        minimizer_matches.write_hits(sequence, *minimizer_hits, k);
                    }
                }

                // infer the clusters of hits
                MinimizerHitClusters clusters_of_hits =
                    get_minimizer_hit_clusters(sample_name, sequence,
                        index.get_prg_max_path_lengths(), index.get_prg_names(),
                        minimizer_hits, max_diff,
                        fraction_kmers_required_for_cluster, cluster_def_file,
                        cluster_filter_file, min_cluster_size,
                        expected_number_kmers_in_read_sketch, rng_seed);

                const std::string sam_record = filtered_mappings.get_sam_record_from_hit_cluster(
                    sequence, clusters_of_hits);

#pragma omp critical(pangraph)
                {
                    add_clusters_to_pangraph(clusters_of_hits, pangraph, index, 0);
                    filtered_mappings.write_sam_record(sam_record);
                    if (!paf_file.is_fake_file) {
                        paf_file.write_clusters(sequence, clusters_of_hits);
                    }
                }
            }

            if (coverageExceeded)
                break; // max_covg exceeded, get out
        }
    }

    BOOST_LOG_TRIVIAL(info) << "Processed " << id << " reads";

    BOOST_LOG_TRIVIAL(debug) << "Pangraph has " << pangraph->nodes.size() << " nodes";

    covg = covg / genome_size;
    BOOST_LOG_TRIVIAL(debug) << "Estimated coverage: " << covg;

    return covg;
}

void open_file_for_reading(const std::string& file_path, std::ifstream& stream)
{
    stream.open(file_path);
    if (!stream.is_open()) {
        fatal_error("Error opening file ", file_path);
    }
}

void open_file_for_writing(const std::string& file_path, std::ofstream& stream)
{
    stream.open(file_path);
    if (!stream.is_open()) {
        fatal_error("Error opening file ", file_path);
    }
}

void open_file_for_appending(const std::string& file_path, std::ofstream& stream)
{
    stream.open(file_path, std::ios::app);
    if (!stream.is_open()) {
        fatal_error("Error opening file ", file_path);
    }
}

uint32_t strtogs(const char* str)
{
    double x;
    char* p;
    x = strtod(str, &p);

    if (x < 0) {
        throw std::logic_error("Negative number passed for genome size");
    }

    if (*p == 'G' || *p == 'g') {
        x *= 1e9;
    } else if (*p == 'M' || *p == 'm') {
        x *= 1e6;
    } else if (*p == 'K' || *p == 'k') {
        x *= 1e3;
    }
    if (x > UINT32_MAX) {
        throw std::runtime_error(
            "Cannot handle genome size larger than 32-bit unsigned integer");
    }
    return (uint32_t)(x + .499);
}

std::string transform_cli_gsize(std::string str)
{
    return int_to_string(strtogs(str.c_str()));
}

std::string make_absolute(std::string str) { return fs::absolute(str).string(); }

std::vector<SampleData> load_read_index(
    const fs::path& read_index_fpath)
{
    std::map<SampleIdText, SampleFpath> samples;
    std::string name, line;
    fs::ifstream instream(read_index_fpath);
    if (instream.fail()) {
        fatal_error("Unable to open read index file ", read_index_fpath);
    }
    while (getline(instream, line)) {
        std::istringstream linestream(line);
        if (std::getline(linestream, name, '\t')) {
            if (samples.find(name) != samples.end()) {
                BOOST_LOG_TRIVIAL(warning)
                    << "Warning: non-unique sample ids given! Only the last "
                       "of these will be kept";
            }
            std::string reads_path;
            linestream >> reads_path;
            if (reads_path.empty()) {
                fatal_error("Malformatted read index file entry for ", name);
            }
            samples[name] = reads_path;
        }
    }
    BOOST_LOG_TRIVIAL(info) << "Finished loading " << samples.size()
                            << " samples from read index";
    return std::vector<SampleData>(
        samples.begin(), samples.end());
}

std::pair<int, std::string> build_memfd(const std::string &data) {
    int fd = memfd_create("pandora_memfd", MFD_ALLOW_SEALING);
    if (fd == -1)
        fatal_error("memfd could not be created");

    /* set the size of the file */
    if (ftruncate(fd, data.length()) == -1)
        fatal_error("Could not truncate memfd");

    if (write(fd, data.c_str(), data.size()) != (long)(data.size()))
        fatal_error("Could not write all the data to memfd");

    if (fcntl(fd, F_ADD_SEALS, F_SEAL_WRITE) == -1)
        fatal_error("Could not add write seal to memfd");

    if (fsync(fd) == -1)
        fatal_error("Could not fsync memfd.");

    std::stringstream ss_filepath;
    ss_filepath << "/proc/" << getpid() << "/fd/" << fd;
    return std::make_pair(fd, ss_filepath.str());
}

void build_file(const std::string &filepath, const std::string &data) {
    std::ofstream output_file;
    open_file_for_writing(filepath, output_file);
    output_file.write(data.c_str(), data.size());
    output_file.close();
}

void concatenate_text_files(
    const fs::path& output_filename, const std::vector<fs::path>& input_filenames,
    const std::string &prepend)
{
    std::ofstream output_filehandler;
    open_file_for_writing(output_filename.string(), output_filehandler);

    if (!prepend.empty()) {
        output_filehandler << prepend << std::endl;
    }

    for (const fs::path& input_filename : input_filenames) {
        std::ifstream input_filehandler;
        open_file_for_reading(input_filename.string(), input_filehandler);
        output_filehandler << input_filehandler.rdbuf();
        input_filehandler.close();
    }

    output_filehandler.close();
}

std::string reverse_complement(const std::string& forward)
{
    const auto len { forward.size() };
    std::string reverse(len, ' ');
    for (size_t k = 0; k < len; k++) {
        const char base { forward[k] };
        const char magic = base & 2 ? 4 : 21;
        reverse[len - k - 1] = base ^ magic;
    }
    reverse[len] = '\0';
    return reverse;
}

std::pair<std::vector<std::string>, std::vector<size_t>> split_ambiguous(const std::string& input_string, uint8_t delim)
{
    std::vector<std::string> substrs;
    std::vector<size_t> offsets;
    auto start { 0 };
    auto current_index { 0 };
    auto valid_substring_length { 0 };
    for (const auto& base : input_string) {
        const uint32_t coded_base = pandora::nt4(base);
        const bool is_ambiguous = coded_base == delim;
        if (is_ambiguous) {
            if (valid_substring_length > 0) {
                substrs.emplace_back(input_string.substr(start, valid_substring_length));
                offsets.emplace_back(start);
            }
            start = current_index + 1;
            valid_substring_length = 0;
        } else {
            ++valid_substring_length;
        }
        ++current_index;
    }
    if (valid_substring_length > 0) {
        substrs.emplace_back(input_string.substr(start, valid_substring_length));
        offsets.emplace_back(start);
    }
    return std::make_pair(substrs, offsets);
}
