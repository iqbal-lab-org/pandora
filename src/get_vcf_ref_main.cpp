#include "get_vcf_ref_main.h"
#include "cli_helpers.h"

void setup_get_vcf_ref_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<GetVcfRefOptions>();

    std::string description
        = "Outputs a fasta suitable for use as the VCF reference using input sequences";
    auto* gvr_subcmd = app.add_subcommand("get_vcf_ref", description);

    gvr_subcmd->add_option("<INDEX>", opt->index_file, "A pandora index (.panidx.zip) file")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->check(PandoraIndexValidator())
        ->transform(make_absolute)
        ->type_name("FILE");

    gvr_subcmd
        ->add_option("<QUERY>", opt->seqfile,
            "Fast{a,q} file of sequences to retrieve the PRG reference for")
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    gvr_subcmd->add_flag(
        "-z,--compress", opt->compress, "Compress the output with gzip");

    gvr_subcmd->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    gvr_subcmd->callback([opt]() { pandora_get_vcf_ref(*opt); });
}

int pandora_get_vcf_ref(GetVcfRefOptions const& opt)
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
    index.load_all_prgs();
    BOOST_LOG_TRIVIAL(info) << "Index loaded successfully!";

    Fastaq output_fasta(opt.compress, false);

    if (opt.seqfile.empty()) {
        for (const auto& prg_ptr : index.get_loaded_prgs()) {
            std::vector<LocalNodePtr> npath;
            npath = prg_ptr->prg.top_path();
            output_fasta.add_entry(prg_ptr->name, prg_ptr->string_along_path(npath));
        }
    } else {
        std::vector<LocalNodePtr> npath;
        std::string read_string;
        FastaqHandler readfile(opt.seqfile.string());
        bool found { false };

        for (const auto& prg_ptr : index.get_loaded_prgs()) {
            found = false;
            readfile.get_nth_read(0);
            while (not readfile.eof()) {
                try {
                    readfile.get_next();
                } catch (std::out_of_range& err) {
                    break;
                }
                npath = prg_ptr->get_valid_vcf_reference(readfile.read);
                if (not npath.empty()) {
                    BOOST_LOG_TRIVIAL(trace) << ">" << prg_ptr->name << std::endl
                                             << readfile.read;
                    output_fasta.add_entry(
                        prg_ptr->name, prg_ptr->string_along_path(npath));
                    found = true;
                    break;
                }
            }

            if (!found) {
                if (npath.empty()) {
                    fatal_error("PRG is empty");
                }
                BOOST_LOG_TRIVIAL(debug)
                    << "Using top path as ref for " << prg_ptr->name;
                npath = prg_ptr->prg.top_path();
                output_fasta.add_entry(
                    prg_ptr->name, prg_ptr->string_along_path(npath));
            }
        }
    }

    std::string prg_file_prefix = opt.index_file.string();
    // removes ".panidx.zip"
    prg_file_prefix = prg_file_prefix.substr(0, prg_file_prefix.size()-11);
    output_fasta.save(prg_file_prefix + ".vcf_ref.fa.gz");

    return 0;
}
