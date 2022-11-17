#include "index_main.h"
#include "localPRG_reader.h"

void setup_index_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<IndexOptions>();
    auto* index_subcmd = app.add_subcommand(
        "index", "Index population reference graph (PRG) sequences.");
    index_subcmd->add_option("<PRG>", opt->prgfile, "PRG to index (in fasta format)")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    index_subcmd
        ->add_option(
            "-w", opt->window_size, "Window size for (w,k)-minimizers (must be <=k)")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd->add_option("-k", opt->kmer_size, "K-mer size for (w,k)-minimizers")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd->add_option("-m", opt->indexing_upper_bound,
                    "Maximum number of minimising kmers or starting walks or nodes to explore per locus "
                    "during indexing. If exceeded, the locus is assessed as too complex and is ignored, which might flag an "
                    "issue with that specific locus. This is used as an upper bound to ignore overly complex PRGs, which "
                    "very likely does not represent the real biological complexity and should be revisited.")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd->add_option("-o,--outfile", opt->outfile, "Filename for the index. Must end in .zip")
        ->type_name("FILE")
        ->check(CLI::NonexistentPath)
        ->check(check_if_is_zip_file)
        ->transform(make_absolute)
        ->default_str("<PRG>.zip");

    index_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    index_subcmd->callback([opt]() { pandora_index(*opt); });
}

int pandora_index(IndexOptions const& opt)
{
    auto log_level = boost::log::trivial::info;
    if (opt.verbosity == 1) {
        log_level = boost::log::trivial::debug;
    } else if (opt.verbosity > 1) {
        log_level = boost::log::trivial::trace;
    }
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= log_level);

    if (opt.window_size > opt.kmer_size) {
        throw std::logic_error("W must NOT be greater than K");
    }
    if (opt.window_size <= 0) {
        throw std::logic_error("W must be a positive integer");
    }
    if (opt.kmer_size <= 0) {
        throw std::logic_error("K must be a positive integer");
    }

    fs::path outfile = opt.outfile;
    if (outfile.empty()) {
        outfile = fs::path(opt.prgfile.string() + ".panidx.zip");
    }

    LocalPRG::do_path_memoization_in_nodes_along_path_method = true;


    BOOST_LOG_TRIVIAL(info) << "Indexing PRG...";
    Index::build_index_on_disk(opt.window_size, opt.kmer_size, opt.prgfile, outfile,
        opt.indexing_upper_bound, opt.threads);
    BOOST_LOG_TRIVIAL(info) << "All done!";

    return 0;
}
