#include "index_main.h"

void setup_index_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<IndexOptions>();
    auto index_subcmd = app.add_subcommand(
        "index", "Index population reference graph (PRG) sequences.");
    index_subcmd->add_option("<PRG>", opt->prgfile, "PRG to index (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    std::stringstream desc;
    desc << "Window size for (w,k)-minimizers (must be <=K) [default: "
         << opt->window_size << "]";
    index_subcmd->add_option("-w", opt->window_size, desc.str())->type_name("INT");

    desc.str(std::string());
    desc << "K-mer size for (w,k)-minimizers [default: " << opt->kmer_size << "]";
    index_subcmd->add_option("-k", opt->kmer_size, desc.str())->type_name("INT");

    desc.str(std::string());
    desc << "Maximum number of threads to use [default: " << opt->threads << "]";
    index_subcmd->add_option("-t,--threads", opt->threads, desc.str())
        ->type_name("INT");

    desc.str(std::string());
    desc << "Offset for PRG ids [default: " << opt->id_offset << "]";
    index_subcmd->add_option("--offset", opt->id_offset, desc.str())->type_name("INT");

    desc.str(std::string());
    desc << "Filename for the index [default: <PRG>.k" << opt->kmer_size << ".w"
         << opt->window_size << ".idx]";
    index_subcmd->add_option("-o,--outfile", opt->outfile, desc.str())
        ->type_name("FILE");

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

    LocalPRG::do_path_memoization_in_nodes_along_path_method = true;

    // load PRGs from file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile, opt.id_offset);

    // get output directory for the gfa
    boost::filesystem::path p(opt.prgfile);
    boost::filesystem::path dir = p.parent_path();
    std::string outdir = dir.string();
    if (outdir.empty())
        outdir = ".";
    outdir += "/kmer_prgs";

    BOOST_LOG_TRIVIAL(info) << "Indexing PRG...";
    auto index = std::make_shared<Index>();
    index_prgs(prgs, index, opt.window_size, opt.kmer_size, outdir, opt.threads);

    // save index
    BOOST_LOG_TRIVIAL(info) << "Saving index...";
    if (not opt.outfile.empty())
        index->save(opt.outfile);
    else if (opt.id_offset > 0)
        index->save(opt.prgfile + "." + std::to_string(opt.id_offset), opt.window_size,
            opt.kmer_size);
    else
        index->save(opt.prgfile, opt.window_size, opt.kmer_size);

    BOOST_LOG_TRIVIAL(info) << "All done!";
    return 0;
}
