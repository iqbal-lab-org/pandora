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
#include "noise_filtering.h"
#include "minihit.h"
#include "fastaq_handler.h"
#include "minimizermatch_file.h"
#include "sam_file.h"

#define _GNU_SOURCE
#include <sys/mman.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <array>
#include <random>

#include <cstdio>
#include <sys/wait.h>

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

void read_prg_file(
    std::vector<std::shared_ptr<LocalPRG>>& prgs, const fs::path& filepath, uint32_t id)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading PRGs from file " << filepath;

    FastaqHandler fh(filepath.string());
    while (!fh.eof()) {
        try {
            fh.get_next();
        } catch (std::out_of_range& err) {
            break;
        }
        if (fh.name.empty() or fh.read.empty())
            continue;
        auto s = std::make_shared<LocalPRG>(LocalPRG(id, fh.name,
            fh.read)); // build a node in the graph, which will represent a LocalPRG
                       // (the graph is a list of nodes, each representing a LocalPRG)
        if (s != nullptr) {
            prgs.push_back(s);
            id++;
        } else {
            fatal_error("Failed to make LocalPRG for ", fh.name);
        }
    }
    BOOST_LOG_TRIVIAL(debug) << "Number of LocalPRGs read: " << prgs.size();
}

void load_PRG_kmergraphs(std::vector<std::shared_ptr<LocalPRG>>& prgs,
    const uint32_t& w, const uint32_t& k, const fs::path& prgfile)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading kmer_prgs from files";
    const auto kmer_prgs_dir { prgfile.parent_path() / "kmer_prgs" };

    auto dir_num = 0;
    fs::path dir;
    for (const auto& prg : prgs) {
        if (prg->id % 4000 == 0) {
            dir = kmer_prgs_dir / int_to_string(dir_num + 1);
            dir_num++;
            if (not fs::exists(dir))
                dir = kmer_prgs_dir;
        }
        const auto filename { prg->name + ".k" + std::to_string(k) + ".w"
            + std::to_string(w) + ".gfa" };
        prg->kmer_prg.load(dir / filename);
    }
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
                minimizer_hits->add_hit(sequence.id, *sequenceSketchIt, miniRecord);
            }
        }
    }
}


void decide_if_add_cluster_or_not(
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs,
    const std::set<MinimizerHitPtr, pComp>::iterator &mh_previous,
    const uint32_t expected_number_kmers_in_read_sketch,
    const float fraction_kmers_required_for_cluster,
    const uint32_t min_cluster_size,
    const Hits &current_cluster,
    const uint32_t number_of_equal_read_minimizers,
    const std::vector<uint32_t> &distances_between_hits,
    ClusterDefFile& cluster_def_file) {
    // keep clusters which cover at least 1/2 the expected number of minihits
    const uint32_t length_based_threshold = std::min(
        prgs[(*mh_previous)->get_prg_id()]->kmer_prg.min_path_length(),
        expected_number_kmers_in_read_sketch)
                                      * fraction_kmers_required_for_cluster;

    const uint32_t number_of_unique_mini_in_cluster = current_cluster.size() - number_of_equal_read_minimizers;
    // const uint32_t cluster_size_threshold = std::max(length_based_threshold, min_cluster_size);
    const uint32_t cluster_size_threshold = min_cluster_size;
    const bool cluster_should_be_accepted = number_of_unique_mini_in_cluster >= cluster_size_threshold;

#pragma omp critical(cluster_def_file)
    {
        cluster_def_file << seq.name << "\t" << prgs[(*mh_previous)->get_prg_id()]->name << "\t";
        if (cluster_should_be_accepted) {
            cluster_def_file << "accepted\t";
            clusters_of_hits.insert(current_cluster);
        }else {
            cluster_def_file << "rejected\t";
        }

        cluster_def_file << current_cluster.size() << "\t"
                         << number_of_equal_read_minimizers << "\t"
                         << number_of_unique_mini_in_cluster << "\t"
                         << length_based_threshold << "\t"
                         << min_cluster_size << "\t";

        for (const auto &distance : distances_between_hits) {
            cluster_def_file << distance << ",";
        }
        cluster_def_file << "\n";
    }
}

void define_clusters(
    const std::string &sample_name,
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs,
    std::shared_ptr<MinimizerHits> minimizer_hits, const int max_diff,
    const float& fraction_kmers_required_for_cluster, const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch,
    ClusterDefFile& cluster_def_file)
{
    const std::string tag = "[Sample: " + sample_name + ", read index: " + to_string(seq.id) + "]: ";

    BOOST_LOG_TRIVIAL(trace) << tag << "Define clusters of hits from the "
                             << minimizer_hits->hits.size() << " hits";

    if (minimizer_hits->hits.empty()) {
        return;
    }

    // A cluster of hits should match same localPRG, each hit not more than max_diff
    // read bases from the last hit (this last bit is to handle repeat genes).
    auto mh_previous = minimizer_hits->hits.begin();
    Hits current_cluster;
    current_cluster.insert(*mh_previous);
    uint32_t number_of_equal_read_minimizers = 0;
    std::vector<uint32_t> distances_between_hits;
    for (auto mh_current = ++minimizer_hits->hits.begin();
         mh_current != minimizer_hits->hits.end(); ++mh_current) {
        const bool read_minimizer_is_the_same =
            (*mh_current)->get_read_start_position() == (*mh_previous)->get_read_start_position();
        number_of_equal_read_minimizers += (int)read_minimizer_is_the_same;

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
            decide_if_add_cluster_or_not(seq, clusters_of_hits, prgs, mh_previous,
            expected_number_kmers_in_read_sketch, fraction_kmers_required_for_cluster,
            min_cluster_size, current_cluster, number_of_equal_read_minimizers,
            distances_between_hits,cluster_def_file);

            // prepare next cluster
            current_cluster.clear();
            number_of_equal_read_minimizers=0;
            distances_between_hits.clear();

        } else {
            if (!read_minimizer_is_the_same) {
                distances_between_hits.push_back(distance_between_hits);
            }
        }
        current_cluster.insert(*mh_current);
        mh_previous = mh_current;
    }
    decide_if_add_cluster_or_not(seq, clusters_of_hits, prgs, mh_previous,
                                 expected_number_kmers_in_read_sketch, fraction_kmers_required_for_cluster,
                                 min_cluster_size, current_cluster, number_of_equal_read_minimizers,
                                 distances_between_hits, cluster_def_file);

    BOOST_LOG_TRIVIAL(trace) << tag << "Found " << clusters_of_hits.size()
                             << " clusters of hits";
}

void filter_clusters(
    const std::string &sample_name,
    const Seq &seq,
    MinimizerHitClusters& clusters_of_hits,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs,
    ClusterFilterFile& cluster_filter_file)
{
    const std::string tag = "[Sample: " + sample_name + ", read index: " + to_string(seq.id) + "]: ";

    // Next order clusters, choose between those that overlap by too much
    BOOST_LOG_TRIVIAL(trace) << tag << "Filter the " << clusters_of_hits.size()
                             << " clusters of hits";
    if (clusters_of_hits.empty()) {
        return;
    }
    // to do this consider pairs of clusters in turn
    auto c_previous = clusters_of_hits.begin();
    for (auto c_current = ++clusters_of_hits.begin();
         c_current != clusters_of_hits.end(); ++c_current) {
        const bool both_clusters_from_same_read =
            ((*(*c_current).begin())->get_read_id() == (*(*c_previous).begin())->get_read_id());
        const bool same_prg = (*(*c_current).begin())->get_prg_id() == (*(*c_previous).begin())->get_prg_id();
        const bool unconsistent_strands = (*(*c_current).begin())->same_strands() != (*(*c_previous).begin())->same_strands();
        const bool current_cluster_is_contained_in_previous =
            (*--(*c_current).end())->get_read_start_position() <=
            (*--(*c_previous).end())->get_read_start_position();
        const bool one_of_the_clusters_should_be_filtered_out =
            both_clusters_from_same_read && ((same_prg && unconsistent_strands) or
                current_cluster_is_contained_in_previous);
        if (one_of_the_clusters_should_be_filtered_out)
        // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters,
        // but could also impose no more than 2k hits in overlap
        {

            if (c_previous->size() >= c_current->size()) {
#pragma omp critical(cluster_filter_file)
                {
                    cluster_filter_file
                        << seq.name << "\t"
                        << prgs[(*c_current->begin())->get_prg_id()]->name << "\t"
                        << c_current->size() << "\t"
                        << "filtered_out\n";
                }
                clusters_of_hits.erase(c_current);
                c_current = c_previous;

            } else {
#pragma omp critical(cluster_filter_file)
                {
                    cluster_filter_file
                        << seq.name << "\t"
                        << prgs[(*c_previous->begin())->get_prg_id()]->name << "\t"
                        << c_previous->size() << "\t"
                        << "filtered_out\n";
                }
                clusters_of_hits.erase(c_previous);
            }
        }
        c_previous = c_current;
    }

    for (const auto &cluster : clusters_of_hits) {
#pragma omp critical(cluster_filter_file)
        {
            cluster_filter_file
                << seq.name << "\t"
                << prgs[(*cluster.begin())->get_prg_id()]->name << "\t"
                << cluster.size() << "\t"
                << "kept\n";
        }

    }

    BOOST_LOG_TRIVIAL(trace) << tag << "Now have " << clusters_of_hits.size()
                             << " clusters of hits";
}

void filter_clusters2(std::set<Hits, clusterComp>& clusters_of_hits,
    const uint32_t& genome_size)
{
    // Sort clusters by size, and filter out those small clusters which are entirely
    // contained in bigger clusters on reads
    BOOST_LOG_TRIVIAL(trace) << "Filter2 the " << clusters_of_hits.size()
                             << " clusters of hits";
    if (clusters_of_hits.empty()) {
        return;
    }

    std::set<Hits, clusterComp_size> clusters_by_size(
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
}

void add_clusters_to_pangraph(
    const MinimizerHitClusters& minimizer_hit_clusters,
    std::shared_ptr<pangenome::Graph> pangraph,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs)
{
    BOOST_LOG_TRIVIAL(trace) << "Add clusters to PanGraph";
    if (minimizer_hit_clusters.empty()) {
        return;
    }

    // to do this consider pairs of clusters in turn
    for (auto cluster : minimizer_hit_clusters) {
        pangraph->add_hits_between_PRG_and_read(prgs[(*cluster.begin())->get_prg_id()],
            (*cluster.begin())->get_read_id(), cluster);
    }
}

MinimizerHitClusters get_minimizer_hit_clusters(
    const std::string &sample_name,
    const Seq &seq,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs,
    std::shared_ptr<MinimizerHits> minimizer_hits,
    std::shared_ptr<pangenome::Graph> pangraph, const int max_diff,
    const uint32_t& genome_size, const float& fraction_kmers_required_for_cluster,
    ClusterDefFile &cluster_def_file,
    ClusterFilterFile &cluster_filter_file,
    const uint32_t min_cluster_size,
    const uint32_t expected_number_kmers_in_read_sketch)
{
    MinimizerHitClusters minimizer_hit_clusters;

    // this step infers the gene order for a read and adds this to the pangraph
    // by defining clusters of hits, keeping those which are not noise and
    // then adding the inferred gene ordering
    if (minimizer_hits->hits.empty()) {
        return minimizer_hit_clusters;
    }

    define_clusters(sample_name, seq, minimizer_hit_clusters, prgs, minimizer_hits, max_diff,
        fraction_kmers_required_for_cluster, min_cluster_size,
        expected_number_kmers_in_read_sketch, cluster_def_file);

    filter_clusters(sample_name, seq, minimizer_hit_clusters, prgs, cluster_filter_file);
    // filter_clusters2(clusters_of_hits, genome_size);

    return minimizer_hit_clusters;
}

// TODO: this should be in a constructor of pangenome::Graph or in a factory class
uint32_t pangraph_from_read_file(const SampleData& sample,
    std::shared_ptr<pangenome::Graph> pangraph, std::shared_ptr<Index> index,
    const std::vector<std::shared_ptr<LocalPRG>>& prgs, const uint32_t w,
    const uint32_t k, const int max_diff, const float& e_rate,
    const fs::path& sample_outdir, const uint32_t min_cluster_size,
    const uint32_t genome_size, const bool illumina, const bool clean,
    const uint32_t max_covg, uint32_t threads, const bool keep_extra_debugging_files)
{
    const SampleIdText sample_name = sample.first;
    const SampleFpath sample_filepath = sample.second;
    const std::string tag = "[Sample " + sample_name + "]: ";

    // constant variables
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

    SAMFile filtered_mappings(sample_outdir / (sample_name + ".filtered.sam"), prgs, k*2);

    MinimizerMatchFile minimizer_matches(sample_outdir / (sample_name + ".minimatches"), prgs, !keep_extra_debugging_files);
    PafFile paf_file(sample_outdir / (sample_name + ".minipaf"), prgs, !keep_extra_debugging_files);
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
                auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits());
                add_read_hits(sequence, minimizer_hits, *index);

                // write unfiltered minimizer hits
#pragma omp critical(minimizer_matches)
                {
                    minimizer_matches.write_hits(sequence, minimizer_hits->hits);
                }

                // infer
                MinimizerHitClusters clusters_of_hits =
                    get_minimizer_hit_clusters(sample_name, sequence, prgs, minimizer_hits, pangraph, max_diff,
                    genome_size, fraction_kmers_required_for_cluster, cluster_def_file,
                    cluster_filter_file, min_cluster_size, expected_number_kmers_in_read_sketch);

#pragma omp critical(pangraph)
                {
                    add_clusters_to_pangraph(clusters_of_hits, pangraph, prgs);
                    filtered_mappings.write_sam_record_from_hit_cluster(
                        sequence, clusters_of_hits);
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

    if (illumina and clean) {
        clean_pangraph_with_debruijn_graph(pangraph, 2, 1, illumina);
        BOOST_LOG_TRIVIAL(debug)
            << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes";
    } else if (clean) {
        clean_pangraph_with_debruijn_graph(pangraph, 3, 1, illumina);
        BOOST_LOG_TRIVIAL(debug)
            << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes";
    }

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

// read all strings in the readsFile file and return them as a vector of strings
std::vector<std::string> get_vector_of_strings_from_file(const std::string& file_path)
{
    std::vector<std::string> lines;
    std::string line;

    std::ifstream in_file;
    open_file_for_reading(file_path, in_file);
    while (getline(in_file, line)) {
        if (line.size() > 0)
            lines.push_back(line);
    }
    in_file.close();

    return lines;
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

std::string remove_spaces_from_string(const std::string& str)
{
    std::string to_return;
    to_return.reserve(str.size());
    for (const char c : str) {
        if (c != '-') {
            to_return += c;
        }
    }
    return to_return;
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

// From https://stackoverflow.com/a/478960/5264075
// Exec a command and returns stdout
std::string exec(const char* cmd) {
    BOOST_LOG_TRIVIAL(debug) << "Running " << cmd;
    std::array<char, 4096> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

void build_file(const std::string &filepath, const std::string &data) {
    std::ofstream output_file;
    open_file_for_writing(filepath, output_file);
    output_file.write(data.c_str(), data.size());
    output_file.close();
}

bool tool_exists(const std::string &command) {
    int ret = std::system(command.c_str());
    return WEXITSTATUS(ret) == 0;
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

int random_int()
{
    static std::random_device dev;
    static std::mt19937 rng(dev());
    static std::uniform_int_distribution<std::mt19937::result_type> dist;
    return dist(rng);
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
