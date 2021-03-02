#include "get_vcf_ref_main.h"

void setup_get_vcf_ref_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<GetVcfRefOptions>();

    std::string description
        = "Outputs a fasta suitable for use as the VCF reference using input sequences";
    auto* gvr_subcmd = app.add_subcommand("get_vcf_ref", description);

    gvr_subcmd->add_option("<PRG>", opt->prgfile, "PRG to index (in fasta format)")
        ->required()
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    gvr_subcmd
        ->add_option("<QUERY>", opt->seqfile,
            "Fast{a,q} file of sequences to retrive the PRG reference for")
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
    std::vector<std::shared_ptr<LocalPRG>> prgs;
    read_prg_file(prgs, opt.prgfile);

    Fastaq output_fasta(opt.compress, false);

    if (opt.seqfile.empty()) {
        for (const auto& prg_ptr : prgs) {
            std::vector<LocalNodePtr> npath;
            npath = prg_ptr->prg.top_path();
            output_fasta.add_entry(prg_ptr->name, prg_ptr->string_along_path(npath));
        }
    } else {
        std::vector<LocalNodePtr> npath;
        std::string read_string;
        FastaqHandler readfile(opt.seqfile);
        bool found { false };

        for (const auto& prg_ptr : prgs) {
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

    std::string prg_file(opt.prgfile);
    output_fasta.save(prg_file + ".vcf_ref.fa.gz");

    return 0;
}
