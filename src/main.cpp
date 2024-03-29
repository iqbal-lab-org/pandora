#include <iostream>
#include "CLI11.hpp"
#include "version.h"
#include "index_main.h"
#include "map_main.h"
#include "compare_main.h"
#include "walk_main.h"
#include "seq2path_main.h"
#include "get_vcf_ref_main.h"
#include "random_main.h"
#include "globals.h"
#include "denovo_discovery/discover_main.h"

class MyFormatter : public CLI::Formatter {
public:
    std::string make_option_opts(const CLI::Option* opt) const override
    {
        std::stringstream out;

        if (opt->get_type_size() != 0) {
            if (!opt->get_type_name().empty())
                out << " " << get_label(opt->get_type_name());
            else if (opt->get_expected_min() > 1)
                out << " x " << opt->get_expected();

            if (opt->get_required())
                out << " " << get_label("[required]");
        }
        if (!opt->get_envname().empty())
            out << " (" << get_label("Env") << ":" << opt->get_envname() << ")";
        if (!opt->get_needs().empty()) {
            out << " " << get_label("Needs") << ":";
            for (const CLI::Option* op : opt->get_needs())
                out << " " << op->get_name();
        }
        if (!opt->get_excludes().empty()) {
            out << " " << get_label("Excludes") << ":";
            for (const CLI::Option* op : opt->get_excludes())
                out << " " << op->get_name();
        }
        return out.str();
    }

    std::string make_option_desc(const CLI::Option* opt) const override
    {
        std::stringstream out;
        out << opt->get_description();
        if (!opt->get_default_str().empty()) {
            out << " [default: " << opt->get_default_str() << "]";
        }
        return out.str();
    }
};

int main(int argc, char* argv[])
{
    CLI::App app { "Pandora: Pan-genome inference and genotyping with long noisy or "
                   "short accurate reads." };
    app.formatter(std::make_shared<MyFormatter>());
    app.add_flag_callback("-V,--version", []() {
        std::cout << "pandora version " << PANDORA_VERSION << std::endl;
        throw(CLI::Success {});
    });
    setup_index_subcommand(app);
    setup_map_subcommand(app);
    setup_compare_subcommand(app);
    setup_discover_subcommand(app);
    setup_walk_subcommand(app);
    setup_seq2path_subcommand(app);
    setup_get_vcf_ref_subcommand(app);
    setup_random_subcommand(app);
    app.require_subcommand();

    for (int i=0; i<argc; i++) {
        PandoraGlobals::command_line += std::string(argv[i]) + " ";
    }

    CLI11_PARSE(app, argc, argv);

    return 0;
}
