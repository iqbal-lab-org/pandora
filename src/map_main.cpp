/*
 * C++ Program to find (w,k)-minimizer hits between a fasta of reads and a set of gene PRGs and output a gfa of genes that there is evidence have been traversed and the order they have been seen.
 */
/*
 * To do:
 * 1. Allow reverse complement strand minimizers
 * 2. Make sure mapped region not got large indels from unsupported intervals
 * 3. Command line argument handling for option of dumped_db screening
 * 4. Change structure so not copying minimzers etc, but use pointers
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <set>
#include <algorithm>
#include <map>
#include <cassert>

#include "utils.h"
#include "localPRG.h"
#include "localgraph.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "index.h"
#include "estimate_parameters.h"
#include "noise_filtering.h"

#include "denovo_discovery/extract_reads.h"
#include "denovo_discovery/denovo_discovery.h"


using std::set;
using std::vector;

namespace fs = boost::filesystem;

static void show_map_usage() {
    std::cerr << "Usage: pandora map -p PRG_FILE -r READ_FILE <option(s)>\n"
              << "Options:\n"
              << "\t-h,--help\t\t\tShow this help message\n"
              << "\t-p,--prg_file PRG_FILE\t\tSpecify a fasta-style prg file\n"
              << "\t-r,--read_file READ_FILE\tSpecify a file of reads in fasta format\n"
              << "\t-o,--outdir OUTDIR\tSpecify directory of output\n"
              << "\t-w W\t\t\t\tWindow size for (w,k)-minimizers, must be <=k, default 14\n"
              << "\t-k K\t\t\t\tK-mer size for (w,k)-minimizers, default 15\n"
              << "\t-m,--max_diff INT\t\tMaximum distance between consecutive hits within a cluster, default 500 (bps)\n"
              << "\t-e,--error_rate FLOAT\t\tEstimated error rate for reads, default 0.11\n"
              << "\t--genome_size\tNUM_BP\tEstimated length of genome, used for coverage estimation\n"
              << "\t--output_kg\t\t\tSave kmer graphs with fwd and rev coverage annotations for found localPRGs\n"
              << "\t--output_vcf\t\t\tSave a vcf file for each found localPRG\n"
              << "\t--vcf_refs REF_FASTA\t\tA fasta file with an entry for each LocalPRG giving reference sequence for\n"
              << "\t\t\t\t\tVCF. Must have a perfect match in the graph and the same name as the graph\n"
              << "\t--output_comparison_paths\tSave a fasta file for a random selection of paths through localPRG\n"
              << "\t--output_covgs\tSave a file of covgs for each localPRG present, one number per base of fasta file\n"
              << "\t--output_mapped_read_fa\tSave a file for each gene containing read parts which overlapped it\n"
              << "\t--illumina\t\t\tData is from illumina rather than nanopore, so is shorter with low error rate\n"
              << "\t--clean\t\t\tAdd a step to clean and detangle the pangraph\n"
              << "\t--bin\t\t\tUse binomial model for kmer coverages, default is negative binomial\n"
              << "\t--max_covg\t\t\tMaximum average coverage from reads to accept\n"
              << "\t--genotype\t\t\tAdd extra step to carefully genotype sites\n"
              << "\t--snps_only\t\t\tWhen genotyping, include only snp sites\n"
              << "\t--discover\t\t\tAdd denovo discovery\n"
              << "\t--denovo_kmer_size\t\t\tKmer size to use for denovo discovery\n"
              << "\t--log_level\t\t\tdebug,[info],warning,error\n"
              << std::endl;
}

int pandora_map(int argc, char *argv[]) {
    // if not enough arguments, print usage
    if (argc < 5) {
        std::cerr << "only provided " << argc << "arguments" << std::endl;
        show_map_usage();
        return 1;
    }

    // otherwise, parse the parameters from the command line
    string prgfile, readfile, outdir = "pandora", vcf_refs_file, log_level="info";
    uint32_t w = 14, k = 15, min_cluster_size = 10, genome_size = 5000000, max_covg = 300,
            min_allele_covg_gt = 0, min_total_covg_gt = 0, min_diff_covg_gt = 0, min_kmer_covg=0; // default parameters
    uint16_t confidence_threshold = 1;
    uint_least16_t denovo_kmer_size{11};
    int max_diff = 250;
    float e_rate = 0.11, min_allele_fraction_covg_gt = 0, genotyping_error_rate=0.01;
    bool output_kg = false, output_vcf = false;
    bool output_comparison_paths = false, output_mapped_read_fa = false;
    bool illumina = false, clean = false;
    bool output_covgs = false, bin = false;
    bool genotype = false, snps_only = false, discover_denovo = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_map_usage();
            return 0;
        } else if ((arg == "-p") || (arg == "--prg_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                prgfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--prg_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-r") || (arg == "--read_file")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                readfile = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--read_file option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-o") || (arg == "--outdir")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                outdir = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--outdir option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-w") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                w = (unsigned) atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-w option requires one argument." << std::endl;
                return 1;
            }
        } else if (arg == "-k") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                k = (unsigned) atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "-k option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-m") || (arg == "--max_diff")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_diff = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_diff option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-c") || (arg == "--min_cluster_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_cluster_size = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_cluster_size option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "-e") || (arg == "--error_rate")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                e_rate = static_cast<float>(atof(
                        argv[++i])); // Increment 'i' so we don't get the argument as the next argv[i].
                if (e_rate < 0.01) {
                    illumina = true;
                }
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--error_rate option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genome_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                genome_size = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--genome_size option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--output_kg")) {
            output_kg = true;
        } else if ((arg == "--output_vcf")) {
            output_vcf = true;
        } else if (arg == "--vcf_refs") {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                vcf_refs_file = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--vcf_refs option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--output_covgs")) {
            output_covgs = true;
        } else if ((arg == "--output_comparison_paths")) {
            output_comparison_paths = true;
        } else if ((arg == "--illumina")) {
            illumina = true;
            if (e_rate > 0.05) {
                e_rate = 0.001;
            }
        } else if ((arg == "--output_mapped_read_fa")) {
            output_mapped_read_fa = true;
        } else if ((arg == "--clean")) {
            clean = true;
        } else if ((arg == "--bin")) {
            bin = true;
        } else if ((arg == "--max_covg")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                max_covg = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--max_covg option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_allele_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_allele_covg_gt = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_allele_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_total_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_total_covg_gt = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_total_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_diff_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_diff_covg_gt = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_diff_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--min_allele_fraction_covg_gt")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                min_allele_fraction_covg_gt = static_cast<float>(atof(argv[++i])); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--min_allele_fraction_covg_gt option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotyping_error_rate")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                genotyping_error_rate = static_cast<float>(atof(argv[++i])); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--genotyping_error_rate option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--confidence_threshold")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                confidence_threshold = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--confidence_threshold option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--denovo_kmer_size")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                denovo_kmer_size = atoi(argv[++i]); // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--denovo_kmer_size option requires one argument." << std::endl;
                return 1;
            }
        } else if ((arg == "--genotype")) {
            genotype = true;
        } else if ((arg == "--snps_only")) {
            snps_only = true;
        } else if ((arg == "--discover")) {
            discover_denovo = true;
        } else if ((arg == "--log_level")) {
            if (i + 1 < argc) { // Make sure we aren't at the end of argv!
                log_level = argv[++i]; // Increment 'i' so we don't get the argument as the next argv[i].
            } else { // Uh-oh, there was no argument to the destination option.
                std::cerr << "--log_level option requires one argument." << std::endl;
                return 1;
            }
        } else {
            cerr << argv[i] << " could not be attributed to any parameter" << endl;
        }
    }

    assert(w <= k);
    assert(not prgfile.empty());
    if (snps_only)
        genotype = true;
    if (genotype)
        output_vcf = true;
    if (illumina and e_rate > 0.1) {
        e_rate = 0.001;
    }
    if (illumina and max_diff > 200) {
        max_diff = 2*k + 1;
    }

    //then run the programme...
    cout << "START: " << now() << endl;
    cout << "\nUsing parameters: " << endl;
    cout << "\tprgfile\t\t" << prgfile << endl;
    cout << "\treadfile\t" << readfile << endl;
    cout << "\toutdir\t" << outdir << endl;
    cout << "\tw\t\t" << w << endl;
    cout << "\tk\t\t" << k << endl;
    cout << "\tmax_diff\t" << max_diff << endl;
    cout << "\terror_rate\t" << e_rate << endl;
    cout << "\toutput_kg\t" << output_kg << endl;
    cout << "\toutput_vcf\t" << output_vcf << endl;
    cout << "\tvcf_refs\t" << vcf_refs_file << endl;
    cout << "\toutput_comparison_paths\t" << output_comparison_paths << endl;
    cout << "\toutput_covgs\t" << output_covgs << endl;
    cout << "\toutput_mapped_read_fa\t" << output_mapped_read_fa << endl;
    cout << "\tillumina\t" << illumina << endl;
    cout << "\tclean\t" << clean << endl;
    cout << "\tbin\t" << bin << endl;
    cout << "\tmax_covg\t" << max_covg << endl;
    cout << "\tgenotype\t" << genotype << endl;
    cout << "\tsnps_only\t" << snps_only << endl;
    cout << "\tdiscover\t" << discover_denovo << endl;
    cout << "\tdenovo_kmer_size\t" << denovo_kmer_size << endl;
    cout << "\tlog_level\t" << log_level << endl << endl;

    auto g_log_level{boost::log::trivial::info};
    if (log_level == "debug")
        g_log_level = boost::log::trivial::debug;
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= g_log_level);

    fs::create_directories(outdir);
    if (output_kg)
        fs::create_directories(outdir + "/kmer_graphs");

    cout << now() << "Loading Index and LocalPRGs from file" << endl;
    auto index = std::make_shared<Index>();
    index->load(prgfile, w, k);
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, prgfile);
    load_PRG_kmergraphs(prgs, w, k, prgfile);

    cout << now() << "Constructing pangenome::Graph from read file (this will take a while)" << endl;
    auto minimizer_hits = std::make_shared<MinimizerHits>(MinimizerHits(100000));
    auto pangraph = std::make_shared<pangenome::Graph>(pangenome::Graph());
    uint32_t covg = pangraph_from_read_file(readfile, minimizer_hits, pangraph, index, prgs, w, k, max_diff, e_rate,
                                            min_cluster_size, genome_size, illumina, clean, max_covg);

    cout << now() << "Finished with index, so clear " << endl;
    index->clear();

    cout << now() << "Finished with minihits, so clear " << endl;
    minimizer_hits->clear();

    if (pangraph->nodes.empty()) {
        cout << "Found non of the LocalPRGs in the reads." << endl;
        cout << "FINISH: " << now() << endl;
        return 0;
    }

    cout << now() << "Writing pangenome::Graph to file " << outdir << "/pandora.pangraph.gfa" << endl;
    write_pangraph_gfa(outdir + "/pandora.pangraph.gfa", pangraph);

    cout << now() << "Update LocalPRGs with hits" << endl;
    uint32_t sample_id = 0;
    pangraph->setup_kmergraphs(prgs);
    pangraph->add_hits_to_kmergraphs(prgs);

    cout << now() << "Estimate parameters for kmer graph model" << endl;
    auto exp_depth_covg = estimate_parameters(pangraph, outdir, k, e_rate, covg, bin, sample_id);
    if (min_kmer_covg == 0)
        min_kmer_covg = exp_depth_covg/10;

    cout << now() << "Find PRG paths and write to files:" << endl;

    Fastaq consensus_fq(true, true);
    VCF master_vcf;

    VCFRefs vcf_refs;
    std::string vcf_ref;
    vector<KmerNodePtr> kmp;
    vector<LocalNodePtr> lmp;
    std::set<std::pair<ReadCoordinate, GeneIntervalInfo>> pangraph_coordinate_pairs;

    if (output_vcf and !vcf_refs_file.empty()) {
        vcf_refs.reserve(prgs.size());
        load_vcf_refs_file(vcf_refs_file, vcf_refs);
    }

    for (auto c = pangraph->nodes.begin(); c != pangraph->nodes.end();) {
        if (output_vcf
            and !vcf_refs_file.empty()
            and vcf_refs.find(prgs[c->second->prg_id]->name) != vcf_refs.end()) {
            vcf_ref = vcf_refs[prgs[c->second->prg_id]->name];
        }

        prgs[c->second->prg_id]->add_consensus_path_to_fastaq(consensus_fq, c->second, kmp, lmp, w, bin, covg);

        if (kmp.empty()) {
            c = pangraph->remove_node(c->second);
            continue;
        }

        if (output_kg) {
            c->second->kmer_prg.save(outdir + "/kmer_graphs/" + c->second->get_name() + ".kg.gfa",
                                     prgs[c->second->prg_id]);
        }

        if (output_vcf) {
            prgs[c->second->prg_id]->add_variants_to_vcf(master_vcf, c->second, vcf_ref, kmp, lmp, min_kmer_covg);
        }

        if (discover_denovo) {
            const auto padding_size{denovo_kmer_size * (uint_least16_t)2};
            denovo_discovery::add_pnode_coordinate_pairs(prgs, pangraph_coordinate_pairs, c->second, lmp, kmp, padding_size);
        }
        ++c;
    }
    consensus_fq.save(outdir + "/pandora.consensus.fq.gz");
    if (output_vcf)
        master_vcf.save(outdir + "/pandora_consensus.vcf", true, true, true, true, true, true, true);

    if (pangraph->nodes.empty()) {
        cout << "All nodes which were found have been removed during cleaning. Is your genome_size accurate?"
             << " Genome size is assumed to be " << genome_size << " and can be updated with --genome_size" << endl
             << "FINISH: " << now() << endl;
        return 0;
    }

    if (genotype) {
        std::vector<uint32_t> exp_depth_covgs = {exp_depth_covg};
        master_vcf.genotype(exp_depth_covgs, genotyping_error_rate, confidence_threshold, min_allele_covg_gt, min_allele_fraction_covg_gt,
                            min_total_covg_gt, min_diff_covg_gt, snps_only);
        if (snps_only)
            master_vcf.save(outdir + "/pandora_genotyped.vcf", true, true, true, true, false, false, false);
        else
            master_vcf.save(outdir + "/pandora_genotyped.vcf", true, true, true, true, true, true, true);
    }

    if (discover_denovo) {
        denovo_discovery::find_candidates(pangraph_coordinate_pairs,
                                          readfile,
                                          fs::path(outdir),
                                          e_rate,
                                          denovo_kmer_size);
    }

    if (output_mapped_read_fa)
        pangraph->save_mapped_read_strings(readfile, outdir);

    pangraph->clear();

    cout << "FINISH: " << now() << endl;
    return 0;
}

