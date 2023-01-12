#include "random_main.h"
#include "cli_helpers.h"

void setup_random_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<RandomOptions>();

    std::string description = "Outputs a fasta of random paths through the PRGs";
    auto* random_subcmd = app.add_subcommand("random", description);

    random_subcmd
        ->add_option("<INDEX>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
        ->type_name("FILE");

    random_subcmd->add_option("-n", opt->num_paths, "Number of paths to output")
        ->capture_default_str()
        ->type_name("INT");

    random_subcmd->add_flag(
        "-z,--compress", opt->compress, "Compress the output with gzip");

    random_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    random_subcmd->callback([opt]() { pandora_random_path(*opt); });
}

int pandora_random_path(RandomOptions const& opt)
{
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

    Fastaq fa(opt.compress, false);

    for (const auto& prg_ptr : index.get_prgs()) {
        std::unordered_set<std::string> random_paths;
        auto skip = 0;
        while (random_paths.size() < opt.num_paths and skip < 10) {
            auto spath = prg_ptr->random_path();
            if (random_paths.find(spath) != random_paths.end()) {
                skip += 1;
            } else {
                random_paths.insert(spath);
            }
        }
        uint32_t i = 0;
        for (const auto& path : random_paths) {
            fa.add_entry(prg_ptr->name + "_" + std::to_string(i), path);
            i++;
        }
    }

    fa.save("random_paths.fa.gz");

    return 0;
}
