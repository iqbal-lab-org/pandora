#include <iostream>
#include <cstring>
#include "CLI11.hpp"
#include "index_main.h"
#include "map_main.h"
#include "compare_main.h"
#include "walk_main.h"
#include "check_kmergraph_main.h"
#include "get_vcf_ref_main.h"
#include "random_path_main.h"

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

int pandora_get_vcf_ref(int argc, char* argv[]);

int pandora_random_path(int argc, char* argv[]);

int pandora_merge_index(int argc, char* argv[]);

static int usage()
{
    std::cerr
        << "         merge_index   allows multiple indexes to be merged (no "
           "compatibility check)\n"
        << std::endl;
    return 1;
}

int main(int argc, char* argv[])
{
    CLI::App app { "Pandora: Pan-genome inference and genotyping with long noisy or "
                   "short accurate reads." };
    app.formatter(std::make_shared<MyFormatter>());
    setup_index_subcommand(app);
    setup_map_subcommand(app);
    setup_compare_subcommand(app);
    setup_walk_subcommand(app);
    setup_check_kmergraph_subcommand(app);
    setup_get_vcf_ref_subcommand(app);
    setup_random_path_subcommand(app);
    app.require_subcommand();

    CLI11_PARSE(app, argc, argv);

    return 0;

    //    else if (strcmp(argv[1], "merge_index") == 0)
    //        ret = pandora_merge_index(argc - 1, argv + 1);
}
