#include "walk_main.h"

void setup_walk_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<WalkOptions>();
    auto* walk_subcmd = app.add_subcommand("walk",
        "Outputs a path through the nodes in a PRG corresponding to the either an "
        "input sequence (if it exists) or the top/bottom path");
    walk_subcmd->add_option("<PRG>", opt->prgfile, "A PRG file (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    auto* input = walk_subcmd
                      ->add_option("-i,--input", opt->seqfile,
                          "Fast{a,q} of sequences to output paths through the PRG for")
                      ->check(CLI::ExistingFile.description(""))
                      ->type_name("FILE");

    auto* top = walk_subcmd->add_flag(
        "-T,--top", opt->top, "Output the top path through each local PRG");
    auto* bottom = walk_subcmd->add_flag(
        "-B,--bottom", opt->bottom, "Output the bottom path through each local PRG");

    input->excludes(top)->excludes(bottom);
    top->excludes(bottom);

    // Set the function that will be called when this subcommand is issued.
    walk_subcmd->callback([opt]() { pandora_walk(*opt); });
}

int pandora_walk(WalkOptions const& opt)
{
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);

    std::vector<LocalNodePtr> npath;

    if (opt.top) {
        for (const auto& prg_ptr : prgs) {
            npath = prg_ptr->prg.top_path();
            std::cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
    } else if (opt.bottom) {
        for (const auto& prg_ptr : prgs) {
            npath = prg_ptr->prg.bottom_path();
            std::cout << prg_ptr->name << "\t";
            for (uint32_t j = 0; j != npath.size(); ++j) {
                std::cout << "->" << npath[j]->id;
            }
            std::cout << std::endl;
        }
    } else if (!opt.seqfile.empty()) {

        // for each read in readfile,  infer node path along sequence
        std::string read_string;
        FastaqHandler readfile(opt.seqfile);
        while (not readfile.eof()) {
            try {
                readfile.get_next();
            } catch (std::out_of_range& err) {
                break;
            }

            for (const auto& prg_ptr : prgs) {
                npath = prg_ptr->prg.nodes_along_string(readfile.read);
                if (not npath.empty()) {
                    std::cout << readfile.name << "\t" << prg_ptr->name << "\t";
                    for (uint32_t j = 0; j != npath.size(); ++j) {
                        std::cout << "->" << npath[j]->id;
                    }
                    std::cout << std::endl;
                }
            }
        }
    } else {
        std::cerr << "One of --top, --bottom or --input must be given" << std::endl;
        exit(1);
    }
    return 0;
}
