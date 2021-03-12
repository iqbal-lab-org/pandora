#include "merge_index_main.h"

void setup_merge_index_subcommand(CLI::App& app)
{

    auto opt = std::make_shared<MergeIndexOptions>();

    std::string description
        = "Allows multiple indices to be merged (no compatibility check)";
    auto* merge_subcmd = app.add_subcommand("merge_index", description);

    merge_subcmd->add_option("<IDX>", opt->indicies, "Indices to merge")
        ->required()
        ->type_name("FILES");

    merge_subcmd->add_option("-o,--outfile", opt->outfile, "Filename for merged index")
        ->type_name("FILE")
        ->capture_default_str();

    merge_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    merge_subcmd->callback([opt]() { pandora_merge_index(*opt); });
}

int pandora_merge_index(MergeIndexOptions const& opt)
{
    auto log_level = boost::log::trivial::info;
    if (opt.verbosity == 1) {
        log_level = boost::log::trivial::debug;
    } else if (opt.verbosity > 1) {
        log_level = boost::log::trivial::trace;
    }
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= log_level);

    // merge indexes
    auto index = std::make_shared<Index>();
    for (const auto& indexfile : opt.indicies) {
        index->load(indexfile);
    }

    index->save(opt.outfile.string());

    return 0;
}
