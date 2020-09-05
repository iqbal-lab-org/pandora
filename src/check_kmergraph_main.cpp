#include "check_kmergraph_main.h"

void setup_check_kmergraph_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<CheckKmerGraphOptions>();
    // todo: improve description
    std::string description = "For each sequence, return the path through the PRG";
    auto check_subcmd = app.add_subcommand("check_kmergraph", description);

    check_subcmd->add_option("<PRG>", opt->prgfile, "PRG to index (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    auto input = check_subcmd
                     ->add_option("-i,--input", opt->seqfile,
                         "Fast{a,q} of sequences to output paths through the PRG for")
                     ->check(CLI::ExistingFile.description(""))
                     ->type_name("FILE");

    check_subcmd
        ->add_option(
            "-w", opt->window_size, "Window size for (w,k)-minimizers (must be <=k)")
        ->type_name("INT")
        ->capture_default_str();

    check_subcmd->add_option("-k", opt->kmer_size, "K-mer size for (w,k)-minimizers")
        ->type_name("INT")
        ->capture_default_str();

    auto top = check_subcmd->add_flag(
        "-T,--top", opt->top, "Output the top path through each local PRG");
    auto bottom = check_subcmd->add_flag(
        "-B,--bottom", opt->bottom, "Output the bottom path through each local PRG");

    // todo: this "flag" doesn't seem to do what it says. i.e. it still outputs the node path?
    auto check = check_subcmd->add_flag(
        "--flag", opt->flag, "output success/fail rather than the node path");

    check_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    input->excludes(top)->excludes(bottom);
    check->needs(input);
    top->excludes(bottom);

    check_subcmd->callback([opt]() { pandora_check_kmergraph(*opt); });
}

int pandora_check_kmergraph(CheckKmerGraphOptions const& opt)
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

    // load prg graphs and kmergraphs, for now assume there is only one PRG in this file
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);
    load_PRG_kmergraphs(prgs, opt.window_size, opt.kmer_size, opt.prgfile);

    if (prgs.empty()) {
        BOOST_LOG_TRIVIAL(error) << "PRG is empty!";
        exit(1);
    }

    if (opt.top) {
        for (uint32_t i = 0; i != prgs.size(); ++i) {
            std::cout << "Top node path along PRG " << prgs[i]->name << ": ";
            auto npath = prgs[i]->prg.top_path();
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    } else if (opt.bottom) {
        for (uint32_t i = 0; i != prgs.size(); ++i) {
            std::cout << "Bottom node path along PRG " << prgs[i]->name << ": ";
            auto npath = prgs[i]->prg.bottom_path();
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
        return 0;
    } else if (!opt.seqfile.empty()) {
        // for each read in readfile,  infer node path along sequence
        std::vector<LocalNodePtr> npath;

        FastaqHandler seq_handle(opt.seqfile);

        while (!seq_handle.eof()) {
            try {
                seq_handle.get_next();
            } catch (std::out_of_range& err) {
                break;
            }

            if (prgs.size() == 1) {
                std::cout << "Node path for read " << seq_handle.num_reads_parsed << " "
                          << seq_handle.name << " along PRG " << prgs[0]->name << ": ";
                npath = prgs[0]->prg.nodes_along_string(seq_handle.read);
                if (npath.empty()) {
                    npath = prgs[0]->prg.nodes_along_string(
                        rev_complement(seq_handle.read));
                }
            } else if (seq_handle.num_reads_parsed < prgs.size()) {
                std::cout << "Node path for read " << seq_handle.num_reads_parsed << " "
                          << seq_handle.name << " along PRG "
                          << prgs[seq_handle.num_reads_parsed]->name << ": ";
                npath = prgs[seq_handle.num_reads_parsed]->prg.nodes_along_string(
                    seq_handle.read);
                if (npath.empty()) {
                    npath = prgs[seq_handle.num_reads_parsed]->prg.nodes_along_string(
                        rev_complement(seq_handle.read));
                }
            } else {
                BOOST_LOG_TRIVIAL(error)
                    << "Different numbers of PRGs and reads, exiting";
                exit(1);
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
        BOOST_LOG_TRIVIAL(error) << "One of --top, --bottom or --input must be given";
        exit(1);
    }

    return 0;
}
