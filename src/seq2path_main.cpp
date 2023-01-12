#include "seq2path_main.h"
#include "cli_helpers.h"

void setup_seq2path_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<Seq2PathOptions>();
    std::string description = "For each sequence, return the path through the PRG";
    auto* seq2path_subcmd = app.add_subcommand("seq2path", description);

    seq2path_subcmd->add_option("<INDEX>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
        ->type_name("FILE");

    auto* input = seq2path_subcmd
                      ->add_option("-i,--input", opt->seqfile,
                          "Fast{a,q} of sequences to output paths through the PRG for")
                      ->check(CLI::ExistingFile.description(""))
                      ->type_name("FILE");

    auto* top = seq2path_subcmd->add_flag(
        "-T,--top", opt->top, "Output the top path through each local PRG");
    auto* bottom = seq2path_subcmd->add_flag(
        "-B,--bottom", opt->bottom, "Output the bottom path through each local PRG");

    // todo: this "flag" doesn't seem to do what it says. i.e. it still outputs the node
    // path?
    auto* check = seq2path_subcmd->add_flag(
        "--flag", opt->flag, "output success/fail rather than the node path");

    seq2path_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    input->excludes(top)->excludes(bottom);
    check->needs(input);
    top->excludes(bottom);

    seq2path_subcmd->callback([opt]() { pandora_seq2path(*opt); });
}

int pandora_seq2path(Seq2PathOptions const& opt)
{
    // can either provide a prgfile with 1 prg and a sequence file (or top/bottom) and
    // return the path through the prg for each sequence OR can provide a prgfile with
    // multiple sequences in it, and a seq file with the same number, with 1-1
    // correspondance and check seq n in prg n.
    auto log_level = boost::log::trivial::info;
    if (opt.verbosity == 1) {
        log_level = boost::log::trivial::debug;
    } else if (opt.verbosity > 1) {
        log_level = boost::log::trivial::trace;
    }
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= log_level);

    BOOST_LOG_TRIVIAL(info) << "Loading Index...";
    Index index = Index::load(opt.index_file.string());
    BOOST_LOG_TRIVIAL(info) << "Index loaded successfully!";

    if (index.get_number_of_prgs() == 0) {
        fatal_error("PRG is empty!");
    }

    if (opt.top) {
        for (const auto &prg : index.get_prgs()) {
            std::cout << "Top node path along PRG " << prg->name << ": ";
            auto npath = prg->prg.top_path();
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    } else if (opt.bottom) {
        for (const auto &prg : index.get_prgs()) {
            std::cout << "Bottom node path along PRG " << prg->name << ": ";
            auto npath = prg->prg.bottom_path();
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    } else if (!opt.seqfile.empty()) {
        // for each read in readfile,  infer node path along sequence
        std::vector<LocalNodePtr> npath;

        FastaqHandler seq_handle(opt.seqfile.string());

        while (!seq_handle.eof()) {
            try {
                seq_handle.get_next();
            } catch (std::out_of_range& err) {
                break;
            }

            if (index.get_number_of_prgs() == 1) {
                const auto &prg = index.get_prgs()[0];
                std::cout << "Node path for read " << seq_handle.num_reads_parsed << " "
                          << seq_handle.name << " along PRG " << prg->name << ": ";
                npath = prg->prg.nodes_along_string(seq_handle.read);
                if (npath.empty()) {
                    npath = prg->prg.nodes_along_string(
                        rev_complement(seq_handle.read));
                }
            } else if (seq_handle.num_reads_parsed < index.get_number_of_prgs()) {
                const auto &prg = index.get_prgs()[seq_handle.num_reads_parsed];
                std::cout << "Node path for read " << seq_handle.num_reads_parsed << " "
                          << seq_handle.name << " along PRG "
                          << prg->name << ": ";
                npath = prg->prg.nodes_along_string(seq_handle.read);
                if (npath.empty()) {
                    npath = prg->prg.nodes_along_string(rev_complement(seq_handle.read));
                }
            } else {
                fatal_error("Different numbers of PRGs and reads");
            }
            if (opt.flag) {
                if (npath.empty() and seq_handle.read.size() < 300) {
                    BOOST_LOG_TRIVIAL(error) << "short fail!";
                } else if (npath.empty() and seq_handle.read.size() >= 300) {
                    BOOST_LOG_TRIVIAL(error) << "long fail!";
                } else {
                    BOOST_LOG_TRIVIAL(debug) << "success!";
                }
            } else {
                for (uint32_t j = 0; j != npath.size(); ++j) {
                    std::cout << "->" << *npath[j];
                }
                std::cout << std::endl;
            }
        }
        seq_handle.close();
    } else {
        fatal_error("One of --top, --bottom or --input must be given");
    }

    return 0;
}
