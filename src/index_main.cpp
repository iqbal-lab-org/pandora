#include "index_main.h"

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

    index_subcmd->add_option("-m", opt->max_nb_minimiser_kmers,
                    "Maximum number of minimising kmers per locus. If exceeded, locus "
                    "is ignored")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use")
        ->type_name("INT")
        ->capture_default_str();

    index_subcmd->add_option("-o,--outfile", opt->outfile, "Filename for the index")
        ->type_name("FILE")
        ->transform(make_absolute)
        ->default_str("<PRG>.kXX.wXX.idx");

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

    LocalPRG::do_path_memoization_in_nodes_along_path_method = true;

    // load PRGs from file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile, opt.id_offset);

    // get output directory for the gfa
    const auto kmer_prgs_outdir { opt.prgfile.parent_path() / "kmer_prgs" };

    BOOST_LOG_TRIVIAL(info) << "Indexing PRG...";
    auto index = std::make_shared<Index>();
    index_prgs(
        prgs, index, opt.window_size, opt.kmer_size, opt.max_nb_minimiser_kmers,
        kmer_prgs_outdir, opt.threads);

    // save index
    BOOST_LOG_TRIVIAL(info) << "Saving index...";
    if (not opt.outfile.empty()) {
        index->save(opt.outfile);
    } else if (opt.id_offset > 0) {
        const fs::path outfile { opt.prgfile.string() + "."
            + std::to_string(opt.id_offset) };
        index->save(outfile, opt.window_size, opt.kmer_size);
    } else {
        index->save(opt.prgfile, opt.window_size, opt.kmer_size);
    }

    index->clear();

    BOOST_LOG_TRIVIAL(info) << "All done!";
    return 0;
}
