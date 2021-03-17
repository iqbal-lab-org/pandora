#include "vcf.h"

void VCF::add_record_core(const VCFRecord& vr)
{
    records.push_back(std::make_shared<VCFRecord>(vr));
    chrom_to_record_interval_tree[vr.get_chrom()].add(
        vr.get_pos(), vr.get_ref_end_pos(), records.back().get());
}

void VCF::add_record(const std::string& chrom, uint32_t position,
    const std::string& ref, const std::string& alt, const std::string& info,
    const std::string& graph_type_info)
{
    VCFRecord vr(this, chrom, position, ref, alt, info, graph_type_info);
    this->add_record(vr);
}

void VCF::add_record(const VCFRecord& vcf_record)
{
    if (find_record_in_records(vcf_record)
        == records.end()) { // TODO: improve this search to log(n) using a map or sth
        add_record_core(vcf_record);
    }
}

VCFRecord& VCF::add_or_update_record_restricted_to_the_given_samples(
    VCFRecord& vr, const std::vector<std::string>& sample_names)
{
    // TODO: refactor this, this function does too much
    const bool record_and_samples_are_consistent
        = vr.sampleIndex_to_sampleInfo.size() == sample_names.size()
        or sample_names.size() == 0;
    if (!record_and_samples_are_consistent) {
        fatal_error("Error updating record to a subset of samples: record and subset "
                    "of samples given are inconsistent");
    }

    auto record_it = find_record_in_records(
        vr); // TODO: improve this search to log(n) using a map or sth
    if (record_it == records.end()) {
        add_record_core(vr);
        record_it = --records.end();
        (*record_it)
            ->reset_sample_infos_to_contain_the_given_number_of_samples(samples.size());
    }

    for (uint32_t i = 0; i < sample_names.size(); ++i) {
        auto& name = sample_names[i];
        ptrdiff_t sample_index
            = get_sample_index(name); // TODO: potentially add new samples to the VCF
        (*record_it)->sampleIndex_to_sampleInfo[sample_index]
            = vr.sampleIndex_to_sampleInfo[i]; // TODO: assumes that sample_names
                                               // indexes == vr's sample names
    }

    return **record_it;
}

void VCF::add_samples(const std::vector<std::string>& sample_names)
{
    for (uint32_t i = 0; i < sample_names.size(); ++i) {
        auto& name = sample_names[i];
        get_sample_index(name);
    }
}

ptrdiff_t VCF::get_sample_index(const std::string& name)
{
    // if this sample has not been added before, add a column for it
    auto sample_it = find(samples.begin(), samples.end(), name);
    if (sample_it == samples.end()) {
        samples.push_back(name);
        for (auto& record_ptr : records) {
            record_ptr->add_new_samples(1);

            const bool record_samples_match_VCF_samples
                = samples.size() == record_ptr->sampleIndex_to_sampleInfo.size();
            if (!record_samples_match_VCF_samples) {
                fatal_error(
                    "Error on adding a sample to VCF record: VCF record samples "
                    "do no match global VCF samples");
            }
        }
        return samples.size() - 1;
    } else {
        return distance(samples.begin(), sample_it);
    }
}

// TODO: this method has two responsabilities, split into two
void VCF::add_a_new_record_discovered_in_a_sample_and_genotype_it(
    const std::string& sample_name, const std::string& chrom, const uint32_t pos,
    const std::string& ref, const std::string& alt)
{
    if (ref == "" and alt == "") {
        return;
    }

    ptrdiff_t sample_index = get_sample_index(sample_name);

    VCFRecord vcf_record(this, chrom, pos, ref, alt);
    VCFRecord* vcf_record_pointer;

    auto vcf_record_iterator = find_record_in_records(
        vcf_record); // TODO: improve this search to log(n) using alt map or sth
    const bool vcf_record_was_found = vcf_record_iterator != records.end();
    if (vcf_record_was_found) {
        (*vcf_record_iterator)
            ->sampleIndex_to_sampleInfo[sample_index]
            .set_gt_from_max_likelihood_path(1);
        vcf_record_pointer = vcf_record_iterator->get();
    } else {
        // either we have the ref allele, an alternative allele for alt too nested site,
        // or alt mistake
        bool vcf_record_was_processed = false;

        const bool sample_genotyped_towards_ref_allele = ref == alt;
        if (sample_genotyped_towards_ref_allele) {
            // TODO: create a method to find records based on chrom, pos and ref only
            for (const auto& record : records) {
                if (record->get_chrom() == chrom and record->get_pos() == pos
                    and record->get_ref() == ref) {
                    record->sampleIndex_to_sampleInfo[sample_index]
                        .set_gt_from_max_likelihood_path(0);
                    vcf_record_pointer = record.get();
                    vcf_record_was_processed = true;
                }
            }
        } else {
            add_record(
                chrom, pos, ref, alt, "SVTYPE=COMPLEX", "GRAPHTYPE=TOO_MANY_ALTS");
            records.back()
                ->sampleIndex_to_sampleInfo[sample_index]
                .set_gt_from_max_likelihood_path(1);
            vcf_record_pointer = records.back().get();
            vcf_record_was_processed = true;
        }

        // check if there was a mistake
        if (!vcf_record_was_processed) {
            fatal_error("Error when adding a new VCF record discovered in a sample");
        }
    }

    update_other_samples_of_this_record(vcf_record_pointer);
}

void VCF::update_other_samples_of_this_record(VCFRecord* reference_record)
{
    // update other samples at this site if they have ref allele at this pos
    for (const auto& other_record : records) {
        const bool both_records_are_on_the_same_site
            = other_record->get_chrom() == reference_record->get_chrom();
        const bool reference_record_start_overlaps_other_record
            = other_record->get_pos() <= reference_record->get_pos()
            and reference_record->get_pos()
                < other_record->get_pos() + other_record->get_ref().length();
        if (both_records_are_on_the_same_site
            and reference_record_start_overlaps_other_record) {
            for (uint32_t sample_index = 0;
                 sample_index != other_record->sampleIndex_to_sampleInfo.size();
                 ++sample_index) {
                if (other_record->sampleIndex_to_sampleInfo[sample_index]
                        .is_gt_from_max_likelihood_path_valid()
                    and other_record->sampleIndex_to_sampleInfo[sample_index]
                            .get_gt_from_max_likelihood_path()
                        == 0) {
                    reference_record->sampleIndex_to_sampleInfo[sample_index]
                        .set_gt_from_max_likelihood_path(0);
                }
            }
        }
    }
}

void VCF::set_sample_gt_to_ref_allele_for_records_in_the_interval(
    const std::string& sample_name, const std::string& chrom, const uint32_t& pos_from,
    const uint32_t& pos_to)
{

    ptrdiff_t sample_index = get_sample_index(sample_name);

    for (auto& record : records) {
        if (record->ref_allele_is_inside_given_interval(chrom, pos_from, pos_to)) {
            record->sampleIndex_to_sampleInfo[sample_index]
                .set_gt_from_max_likelihood_path(0);
        }
    }
}

void VCF::append_vcf(const VCF& other_vcf)
{
    auto original_size = records.size();
    auto num_samples_added = 0;

    BOOST_LOG_TRIVIAL(debug) << "find which samples are new of the "
                             << other_vcf.samples.size() << " samples";
    std::vector<uint_least16_t> other_sample_positions;
    for (const auto& sample : other_vcf.samples) {
        auto sample_it = find(samples.begin(), samples.end(), sample);
        if (sample_it == samples.end()) {
            samples.push_back(sample);
            other_sample_positions.push_back(samples.size() - 1);
            num_samples_added += 1;
        } else {
            other_sample_positions.push_back(distance(samples.begin(), sample_it));
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "for all existing " << original_size
                             << " records, add null entries for the "
                             << num_samples_added << " new samples";
    for (uint_least64_t i = 0; i < original_size; ++i) {
        records[i]->add_new_samples(num_samples_added);
    }

    BOOST_LOG_TRIVIAL(debug) << "add the " << other_vcf.records.size() << " records";
    for (auto& recordPointer : other_vcf.records) {
        auto& record = *recordPointer;
        VCFRecord& vr = add_or_update_record_restricted_to_the_given_samples(
            record, other_vcf.samples);
        for (uint_least16_t j = 0; j < other_vcf.samples.size(); ++j) {
            vr.sampleIndex_to_sampleInfo[other_sample_positions[j]]
                = record.sampleIndex_to_sampleInfo[j];
            // NB this overwrites old data without checking
        }
    }
}

void VCF::sort_records()
{
    sort(records.begin(), records.end(),
        [](const std::shared_ptr<VCFRecord>& lhs,
            const std::shared_ptr<VCFRecord>& rhs) { return (*lhs) < (*rhs); });
}

bool VCF::pos_in_range(
    const uint32_t from, const uint32_t to, const std::string& chrom) const
{
    // TODO : performance improvement: use this->chrom_to_record_interval_tree

    // is there a record contained in the range from,to?
    for (const auto& recordPointer : records) {
        const auto& record = *recordPointer;
        if (chrom == record.get_chrom() and from < record.get_pos()
            and record.get_pos() + record.get_ref().length() <= to) {
            return true;
        }
    }
    return false;
}

void VCF::genotype(const bool do_local_genotyping)
{
    const bool all_SV_types = not genotyping_options->is_snps_only();

    for (auto& vcf_record : records) {
        const bool should_genotype_record = all_SV_types
            or (genotyping_options->is_snps_only() and vcf_record->is_SNP());
        if (should_genotype_record) {
            if (do_local_genotyping) {
                vcf_record->genotype_from_coverage();
            } else {
                vcf_record
                    ->genotype_from_coverage_using_maximum_likelihood_path_as_reference();
            }
        }
    }

    if (do_local_genotyping) {
        make_gt_compatible();
    }
}

void VCF::merge_multi_allelic_core(VCF& merged_VCF, uint32_t max_allele_length) const
{
    VCF empty_vcf = VCF(merged_VCF.genotyping_options);
    const bool merged_VCF_passed_as_parameter_is_initially_empty
        = merged_VCF == empty_vcf;
    if (!merged_VCF_passed_as_parameter_is_initially_empty) {
        fatal_error("Error on merging VCFs: initial VCF is not empty");
    }

    size_t vcf_size = this->get_VCF_size();
    const bool no_need_for_merging = vcf_size <= 1;
    if (no_need_for_merging) {
        merged_VCF = *this;
        return;
    }

    merged_VCF.add_samples(this->samples);

    std::shared_ptr<VCFRecord> vcf_record_merged
        = (records[0])->make_copy_as_shared_ptr();
    for_each(records.begin() + 1, records.end(),
        [&](const std::shared_ptr<VCFRecord>& vcf_record_to_be_merged_in_pointer) {
            const bool vcf_record_should_be_merged_in
                = vcf_record_merged->can_biallelic_record_be_merged_into_this(
                    *vcf_record_to_be_merged_in_pointer, max_allele_length);

            if (vcf_record_should_be_merged_in) {
                vcf_record_merged->merge_record_into_this(
                    *vcf_record_to_be_merged_in_pointer);
            } else {
                merged_VCF.add_record_core(*vcf_record_merged);
                vcf_record_merged
                    = vcf_record_to_be_merged_in_pointer->make_copy_as_shared_ptr();
            }
        });
    merged_VCF.add_record_core(*vcf_record_merged);

    merged_VCF.sort_records();

    const bool merging_did_not_create_any_record
        = merged_VCF.get_VCF_size() <= vcf_size;
    if (!merging_did_not_create_any_record) {
        fatal_error("Error on merging VCFs: new VCF records were created, whereas "
                    "this should not be the case");
    }
}

VCF VCF::correct_dot_alleles(const std::string& vcf_ref, const std::string& chrom) const
{
    // TODO: optimize this to several loci

    // NB need to merge multiallelic before
    // NB cannot add covgs after

    VCF vcf_with_dot_alleles_corrected(genotyping_options);
    vcf_with_dot_alleles_corrected.add_samples(samples);

    for (auto& recordPointer : records) {
        auto& record = *recordPointer;

        const bool we_are_not_in_the_given_chrom = record.get_chrom() != chrom;

        if (we_are_not_in_the_given_chrom) {
            vcf_with_dot_alleles_corrected.add_record(record);
            continue;
        }

        const bool record_pos_refers_to_an_existing_pos_in_vcf_ref
            = vcf_ref.length() >= record.get_pos();
        if (!record_pos_refers_to_an_existing_pos_in_vcf_ref) {
            fatal_error("When correcting dot alleles, a VCF record has an inexistent "
                        "position (",
                record.get_pos(), ") in VCF ref with length ", vcf_ref.length());
        }
        const bool record_contains_dot_allele = record.contains_dot_allele();
        const bool there_is_a_previous_letter = record.get_pos() > 0;
        const bool there_is_a_next_letter
            = record.get_pos() + record.get_ref().length() + 1 < vcf_ref.length();
        bool record_did_not_contain_dot_allele_or_was_corrected = true;
        if (record_contains_dot_allele and there_is_a_previous_letter) {
            char prev_letter = vcf_ref[record.get_pos() - 1];
            record.correct_dot_alleles_adding_nucleotide_before(prev_letter);
        } else if (record_contains_dot_allele and there_is_a_next_letter) {
            char next_letter;
            if (record.allele_is_dot(record.get_ref())) {
                next_letter = vcf_ref[record.get_pos()];
            } else {
                next_letter = vcf_ref[record.get_pos() + record.get_ref().length()];
            }
            record.correct_dot_alleles_adding_nucleotide_after(next_letter);
        } else if (record_contains_dot_allele) {
            record_did_not_contain_dot_allele_or_was_corrected = false;
        }

        if (record_did_not_contain_dot_allele_or_was_corrected) {
            vcf_with_dot_alleles_corrected.add_record(record);
        }
    }

    vcf_with_dot_alleles_corrected.sort_records();

    const bool correcting_dot_alleles_did_not_create_any_record
        = vcf_with_dot_alleles_corrected.get_VCF_size() <= this->get_VCF_size();
    if (!correcting_dot_alleles_did_not_create_any_record) {
        fatal_error(
            "Error on correcting dot alleles: new VCF records were created, whereas "
            "this should not be the case");
    }

    return vcf_with_dot_alleles_corrected;
}

void VCF::make_gt_compatible()
{
    BOOST_LOG_TRIVIAL(info) << now() << "Make all genotypes compatible";

    for (auto& recordPointer : records) {
        VCFRecord& record = *recordPointer;

        std::vector<VCFRecord*> overlapping_records
            = get_all_records_overlapping_the_given_record(record);

        for (VCFRecord* overlapping_record_ptr : overlapping_records) {
            VCFRecord& overlapping_record = *overlapping_record_ptr;

            bool
                record_starts_at_the_same_position_but_ref_is_smaller_than_overlapping_record_ref
                = record.get_pos() == overlapping_record.get_pos()
                and record.get_ref() < overlapping_record.get_ref();

            const bool this_record_starts_before_the_overlapping_record
                = record.get_pos() < overlapping_record.get_pos()
                or record_starts_at_the_same_position_but_ref_is_smaller_than_overlapping_record_ref;

            if (this_record_starts_before_the_overlapping_record) {
                record.solve_incompatible_gt_conflict_with(overlapping_record);
            }
            // else it was already processed, no need to do it twice
        }
    }
}

std::vector<VCFRecord*> VCF::get_all_records_overlapping_the_given_record(
    const VCFRecord& vcf_record)
{
    IITree<uint32_t, VCFRecord*>& record_interval_tree_for_this_record
        = chrom_to_record_interval_tree[vcf_record.get_chrom()];
    record_interval_tree_for_this_record.index();

    std::vector<size_t> overlaps;
    record_interval_tree_for_this_record.overlap(vcf_record.get_pos(),
        vcf_record.get_pos() + vcf_record.get_ref().length(), overlaps);

    std::vector<VCFRecord*> overlapping_records;
    for (size_t i = 0; i < overlaps.size(); ++i) {
        overlapping_records.push_back(
            record_interval_tree_for_this_record.data(overlaps[i]));
    }

    std::sort(overlapping_records.begin(), overlapping_records.end(),
        [](VCFRecord* lhs, VCFRecord* rhs) { return *lhs < *rhs; });

    return overlapping_records;
}

void VCF::save(const fs::path& filepath, bool genotyping_from_maximum_likelihood,
    bool genotyping_from_coverage, bool output_dot_allele, bool graph_is_simple,
    bool graph_is_nested, bool graph_has_too_many_alts, bool sv_type_is_snp,
    bool sv_type_is_indel, bool sv_type_is_ph_snps, bool sv_type_is_complex)
{
    BOOST_LOG_TRIVIAL(debug) << "Saving VCF to " << filepath;
    fs::ofstream handle;
    handle.open(filepath);
    handle << this->to_string(genotyping_from_maximum_likelihood,
        genotyping_from_coverage, output_dot_allele, graph_is_simple, graph_is_nested,
        graph_has_too_many_alts, sv_type_is_snp, sv_type_is_indel, sv_type_is_ph_snps,
        sv_type_is_complex);
    handle.close();
    BOOST_LOG_TRIVIAL(debug) << "Finished saving " << this->records.size()
                             << " entries to file";
}

std::string VCF::get_current_date() const
{
    time_t t = time(0);
    char mbstr[16];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));
    return std::string(mbstr);
}

std::string VCF::header() const
{
    std::string date = get_current_date();

    std::set<std::string> chroms;
    for (const std::shared_ptr<VCFRecord>& record : records) {
        chroms.insert(record->get_chrom());
    }

    std::string header;
    header.reserve(10000);
    header += "##fileformat=VCFv4.3\n";
    header += "##fileDate==";
    header += date;
    header += "\n##ALT=<ID=SNP,Description=\"SNP\">\n";
    header += "##ALT=<ID=PH_SNPs,Description=\"Phased SNPs\">\n";
    header += "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n";
    header += "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and "
              "indels\">\n";
    header
        += "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n";
    header += "##ALT=<ID=SIMPLE,Description=\"Graph bubble is simple\">\n";
    header += "##ALT=<ID=NESTED,Description=\"Variation site was a nested feature in "
              "the graph\">\n";
    header += "##ALT=<ID=TOO_MANY_ALTS,Description=\"Variation site was a multinested "
              "feature with too many alts to include all in the VCF\">\n";
    header += "##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description=\"Type of graph "
              "feature\">\n";
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    header += "##FORMAT=<ID=MEAN_FWD_COVG,Number=R,Type=Integer,Description=\"Mean "
              "forward coverage\">\n";
    header += "##FORMAT=<ID=MEAN_REV_COVG,Number=R,Type=Integer,Description=\"Mean "
              "reverse coverage\">\n";
    header += "##FORMAT=<ID=MED_FWD_COVG,Number=R,Type=Integer,Description=\"Med "
              "forward coverage\">\n";
    header += "##FORMAT=<ID=MED_REV_COVG,Number=R,Type=Integer,Description=\"Med "
              "reverse coverage\">\n";
    header += "##FORMAT=<ID=SUM_FWD_COVG,Number=R,Type=Integer,Description=\"Sum "
              "forward coverage\">\n";
    header += "##FORMAT=<ID=SUM_REV_COVG,Number=R,Type=Integer,Description=\"Sum "
              "reverse coverage\">\n";
    header += "##FORMAT=<ID=GAPS,Number=R,Type=Float,Description=\"Number of gap "
              "bases\">\n";
    header
        += "##FORMAT=<ID=LIKELIHOOD,Number=R,Type=Float,Description=\"Likelihood\">\n";
    header += "##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description=\"Genotype "
              "confidence\">\n";
    for (const std::string& chrom : chroms) {
        header += "##contig=<ID=" + chrom + ">\n";
    }
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (const std::string& sample : samples) {
        header += "\t" + sample;
    }
    header += "\n";
    return header;
}

std::string VCF::to_string(bool genotyping_from_maximum_likelihood,
    bool genotyping_from_coverage, bool output_dot_allele, bool graph_is_simple,
    bool graph_is_nested, bool graph_has_too_many_alts, bool sv_type_is_snp,
    bool sv_type_is_indel, bool sv_type_is_ph_snps, bool sv_type_is_complex)
{
    const bool only_one_flag_is_set
        = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage))
        == 1;
    if (!only_one_flag_is_set) {
        fatal_error(
            "Error on stringifying VCF record: incompatible genotyping options");
    }

    std::stringstream out;
    out << header();

    // TODO: a side-effect of saving a VCF is sorting it, this might not be desirable
    // TODO: remove this side effect or always keep the VCF sorted
    sort_records();

    for (const auto& record : this->records) {
        const bool record_has_dot_allele_and_should_be_output
            = output_dot_allele and record->contains_dot_allele();

        const bool graph_type_condition_is_satisfied
            = (graph_is_simple and record->graph_type_is_simple())
            or (graph_is_nested and record->graph_type_is_nested())
            or (graph_has_too_many_alts and record->graph_type_has_too_many_alts());
        const bool sv_type_condition_is_satisfied
            = (sv_type_is_snp and record->svtype_is_SNP())
            or (sv_type_is_indel and record->svtype_is_indel())
            or (sv_type_is_ph_snps and record->svtype_is_PH_SNPs())
            or (sv_type_is_complex and record->svtype_is_complex());
        const bool graph_and_sv_type_conditions_are_satisfied
            = graph_type_condition_is_satisfied and sv_type_condition_is_satisfied;

        const bool record_should_be_output = record_has_dot_allele_and_should_be_output
            or graph_and_sv_type_conditions_are_satisfied;

        if (record_should_be_output) {
            out << record->to_string(
                genotyping_from_maximum_likelihood, genotyping_from_coverage)
                << std::endl;
        }
    }

    return out.str();
}

bool VCF::operator==(const VCF& y) const
{
    if (records.size() != y.records.size()) {
        return false;
    }
    for (uint32_t i = 0; i != y.records.size(); ++i) {
        if (find_record_in_records(*(y.records[i])) == records.end()) {
            return false;
        }
    }
    return true;
}

bool VCF::operator!=(const VCF& y) const { return !(*this == y); }

void VCF::concatenate_VCFs(const std::vector<fs::path>& VCF_paths_to_be_concatenated,
    const fs::path& final_VCF_file)
{
    std::ofstream out_file;
    open_file_for_writing(final_VCF_file.string(), out_file);

    bool header_is_output = false;
    for (const auto& VCF_path : VCF_paths_to_be_concatenated) {
        std::ifstream in_file;
        open_file_for_reading(VCF_path.string(), in_file);

        std::string line;
        while (std::getline(in_file, line)) {
            if ((line[0] != '#') || (line[0] == '#' && !header_is_output))
                out_file << line << std::endl;
        }
        header_is_output = true;

        in_file.close();
    }
    out_file.close();
}

// find a VCRRecord in records
std::vector<std::shared_ptr<VCFRecord>>::iterator VCF::find_record_in_records(
    const VCFRecord& vr)
{
    return find_if(records.begin(), records.end(),
        [&vr](const std::shared_ptr<VCFRecord>& record) { return *record == vr; });
}
std::vector<std::shared_ptr<VCFRecord>>::const_iterator VCF::find_record_in_records(
    const VCFRecord& vr) const
{
    return find_if(records.begin(), records.end(),
        [&vr](const std::shared_ptr<VCFRecord>& record) { return *record == vr; });
}
