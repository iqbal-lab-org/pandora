#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <numeric>
#include <vector>
#include <algorithm>

#include <boost/log/trivial.hpp>

#include "vcfrecord.h"
#include "vcf.h"
#include "utils.h"
#include "localnode.h"


#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)


void VCF::add_record_core(const VCFRecord &vr) {
    records.push_back(std::make_shared<VCFRecord>(vr));
    chrom2recordIntervalTree[vr.chrom].add(vr.pos, vr.pos + vr.ref.length() + 1, records.back().get());
}

void VCF::add_record(std::string c, uint32_t p, std::string r, std::string a, std::string i, std::string g) {
    VCFRecord vr(c, p, r, a, i, g);
    if (find_record_in_records(vr) == records.end()) { //TODO: improve this search to log(n) using a map or sth
        add_record_core(vr);
        records.back()->sampleIndex_to_sampleInfo.push_back_several_empty_sample_infos(samples.size(), genotyping_options);
    }
}



VCFRecord &VCF::add_record(VCFRecord &vr, const std::vector<std::string> &sample_names) {
    assert(vr.sampleIndex_to_sampleInfo.size() == sample_names.size() or sample_names.size() == 0);

    auto record_it = find_record_in_records(vr); //TODO: improve this search to log(n) using a map or sth
    if (record_it == records.end()) {
        add_record_core(vr);
        records.back()->sampleIndex_to_sampleInfo.clear();
        records.back()->sampleIndex_to_sampleInfo.push_back_several_empty_sample_infos(samples.size(), genotyping_options);
        record_it = --records.end();
    }

    for (uint32_t i=0; i < sample_names.size(); ++i){
        auto &name = sample_names[i];
        ptrdiff_t sample_index = get_sample_index(name);
        (*record_it)->sampleIndex_to_sampleInfo[sample_index] = vr.sampleIndex_to_sampleInfo[i];
    }

    return **record_it;
}

void VCF::add_samples(const std::vector<std::string> sample_names) {
    for (uint32_t i=0; i < sample_names.size(); ++i){
        auto &name = sample_names[i];
        ptrdiff_t sample_index = get_sample_index(name);
    }
}

ptrdiff_t VCF::get_sample_index(const std::string &name) {
    // if this sample has not been added before, add a column for it
    auto sample_it = find(samples.begin(), samples.end(), name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(name);
        for (uint32_t i = 0; i != records.size(); ++i) {
            records[i]->sampleIndex_to_sampleInfo.push_back_several_empty_sample_infos(1, genotyping_options);
            assert(samples.size() == records[i]->sampleIndex_to_sampleInfo.size());
        }
        return samples.size() - 1;
    } else {
        return distance(samples.begin(), sample_it);
    }
}

// TODO: this method has two responsabilities, split into two
void VCF::add_a_new_record_discovered_in_a_sample_and_genotype_it(const std::string &sample_name, const std::string &chrom, const uint32_t pos, const std::string &ref,
                                                                  const std::string &alt) {
    //cout << "adding gt " << chrom << " " << pos << " " << ref << " vs " << sample_name << " " << alt << endl;
    if (ref == "" and alt == "") {
        return;
    }

    //cout << "adding gt " << ref << " vs " << alt << endl;

    ptrdiff_t sample_index = get_sample_index(sample_name);

    VCFRecord vcf_record(chrom, pos, ref, alt);
    VCFRecord *vcf_record_pointer;

    auto vcf_record_iterator = find_record_in_records(vcf_record); //TODO: improve this search to log(n) using alt map or sth
    bool vcf_record_was_found = vcf_record_iterator != records.end();
    if (vcf_record_was_found) {
        //cout << "found record with this ref and alt" << endl;
        (*vcf_record_iterator)->sampleIndex_to_sampleInfo[sample_index].set_gt_from_max_likelihood_path(1);
        vcf_record_pointer = vcf_record_iterator->get();
    } else {
        //cout << "didn't find alt record for pos " << pos << " ref " << ref << " and alt " << alt << endl;
        // either we have the ref allele, an alternative allele for alt too nested site, or alt mistake
        bool vcf_record_was_processed = false;

        bool sample_genotyped_towards_ref_allele = ref == alt;
        if (sample_genotyped_towards_ref_allele) {
            // TODO: create a method to find records based on chrom, pos and ref only
            for (const auto &record : records) {
                if (record->chrom == chrom and record->pos == pos and record->ref == ref) {
                    record->sampleIndex_to_sampleInfo[sample_index].set_gt_from_max_likelihood_path(0);
                    vcf_record_pointer = record.get();
                    vcf_record_was_processed = true;
                }
            }
        }else {
            //cout << "have new allele not in graph" << endl;
            add_record(chrom, pos, ref, alt, "SVTYPE=COMPLEX", "GRAPHTYPE=TOO_MANY_ALTS");
            records.back()->sampleIndex_to_sampleInfo[sample_index].set_gt_from_max_likelihood_path(1);
            vcf_record_pointer = records.back().get();
            vcf_record_was_processed = true;
        }

        // check not mistake
        assert(vcf_record_was_processed);
    }

    update_other_samples_of_this_record(vcf_record_pointer);
}

void VCF::update_other_samples_of_this_record(VCFRecord *reference_record) {
    // update other samples at this site if they have ref allele at this pos
    for (const auto &other_record : records) {
        bool both_records_are_on_the_same_site = other_record->chrom == reference_record->chrom;
        bool reference_record_start_overlaps_other_record = other_record->pos <= reference_record->pos and reference_record->pos < other_record->pos + other_record->ref.length();
        if ( both_records_are_on_the_same_site and reference_record_start_overlaps_other_record) {
            for (uint32_t sample_index = 0; sample_index != other_record->sampleIndex_to_sampleInfo.size(); ++sample_index) {
                if (other_record->sampleIndex_to_sampleInfo[sample_index].is_gt_from_max_likelihood_path_valid()
                    and other_record->sampleIndex_to_sampleInfo[sample_index].get_gt_from_max_likelihood_path() == 0) {
                    //cout << "update my record to have ref allele also for sample " << sample_index << endl;
                    //cout << records[record_index] << endl;
                    reference_record->sampleIndex_to_sampleInfo[sample_index].set_gt_from_max_likelihood_path(0);
                }
            }
        }
    }
}


void VCF::set_sample_gt_to_ref_allele_for_records_in_the_interval(const std::string &sample_name, const std::string &chrom, const uint32_t &pos_from,
                                                                  const uint32_t &pos_to) {

    ptrdiff_t sample_index = get_sample_index(sample_name);

    for (auto &record : records) {
        if (record->ref_allele_is_inside_given_interval(chrom, pos_from, pos_to)) {
            record->sampleIndex_to_sampleInfo[sample_index].set_gt_from_max_likelihood_path(0);
            //cout << "update record " << records[i] << endl;
        }
    }
}

void VCF::append_vcf(const VCF &other_vcf) {
    auto original_size = records.size();
    auto num_samples_added = 0;

    BOOST_LOG_TRIVIAL(debug) << "find which samples are new of the " << other_vcf.samples.size() << " samples";
    std::vector<uint_least16_t> other_sample_positions;
    for (const auto &sample : other_vcf.samples) {
        auto sample_it = find(samples.begin(), samples.end(), sample);
        if (sample_it == samples.end()) {
            samples.push_back(sample);
            other_sample_positions.push_back(samples.size() - 1);
            num_samples_added += 1;
        } else {
            other_sample_positions.push_back(distance(samples.begin(), sample_it));
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "for all existing " << original_size << " records, add null entries for the "
                             << num_samples_added << " new samples";
    assert(original_size < std::numeric_limits<uint_least64_t>::max() ||
           assert_msg("VCF size has got too big to use the append feature"));
    for (uint_least64_t i = 0; i < original_size; ++i) {
        records[i]->sampleIndex_to_sampleInfo.push_back_several_empty_sample_infos(num_samples_added, genotyping_options);
    }

    BOOST_LOG_TRIVIAL(debug) << "add the " << other_vcf.records.size() << " records";
    for (auto &recordPointer : other_vcf.records) {
        auto &record = *recordPointer;
        VCFRecord &vr = add_record(record, other_vcf.samples);
        for (uint_least16_t j = 0; j < other_vcf.samples.size(); ++j) {
            vr.sampleIndex_to_sampleInfo[other_sample_positions[j]] = record.sampleIndex_to_sampleInfo[j];
            //NB this overwrites old data without checking
        }
    }
}

void VCF::sort_records() {
    sort(records.begin(), records.end(), [](const std::shared_ptr<VCFRecord>& lhs, const std::shared_ptr<VCFRecord>& rhs) {
        return (*lhs) < (*rhs);
    });
}

bool VCF::pos_in_range(const uint32_t from, const uint32_t to, const std::string &chrom) const {
    // is there a record contained in the range from,to?
    for (const auto &recordPointer : records) {
        const auto &record = *recordPointer;
        if (chrom == record.chrom and from < record.pos and record.pos + record.ref.length() <= to) {
            return true;
        }
    }
    return false;
}

void VCF::genotype() {
    BOOST_LOG_TRIVIAL(info) << now() << "Genotype VCF";
    bool all_SV_types = not genotyping_options->is_snps_only();

    for (auto &vcf_record : records) {
        if (all_SV_types or (genotyping_options->is_snps_only() and vcf_record->is_SNP())) {
            vcf_record->genotype();
        }
    }
    make_gt_compatible();
}

void VCF::clean() {
    VCFRecord dummy;
    for (auto record_it = records.begin(); record_it != records.end();) {
        if (**record_it == dummy)
            record_it = records.erase(record_it);
        else
            record_it++;
    }
}




VCF VCF::merge_multi_allelic(uint32_t max_allele_length) const {
    size_t vcf_size = this->get_VCF_size();
    bool no_need_for_merging = vcf_size <= 1;
    if (no_need_for_merging)
        return VCF(*this);

    VCF merged_VCF(genotyping_options);
    merged_VCF.add_samples(this->samples);

    VCFRecord vcf_record_merged(**(records.begin()));
    for_each(records.begin()+1, records.end(), [&](const std::shared_ptr<VCFRecord> &vcf_record_to_be_merged_in_pointer) {
        const VCFRecord &vcf_record_to_be_merged_in = *vcf_record_to_be_merged_in_pointer;

        bool ensure_we_are_merging_only_biallelic_records = vcf_record_to_be_merged_in.alts.size() == 1;
        assert(ensure_we_are_merging_only_biallelic_records);

        bool both_records_have_the_same_ref = vcf_record_merged.has_non_null_reference()
                                              and vcf_record_to_be_merged_in.has_non_null_reference()
                                              and vcf_record_merged.ref == vcf_record_to_be_merged_in.ref;

        bool all_alleles_have_at_most_max_allele_length =
                vcf_record_merged.get_longest_allele_length() <= max_allele_length
                and vcf_record_to_be_merged_in.get_longest_allele_length() <= max_allele_length;

        bool vcf_record_should_be_merged_in = vcf_record_to_be_merged_in != vcf_record_merged
                                              and vcf_record_merged.has_the_same_position(vcf_record_to_be_merged_in)
                                              and both_records_have_the_same_ref
                                              and all_alleles_have_at_most_max_allele_length;

        if (vcf_record_should_be_merged_in) {
            // TODO: this code is not covered by tests, IDK what it is supposed to do - commenting it out and asserting out if we reach it
            if (vcf_record_to_be_merged_in.sampleIndex_to_sampleInfo.empty()) {
                assert_msg("VCF::merge_multi_allelic: vcf_record_to_be_merged_in has no samples");
                /*
                records.at(current_vcf_record_position)->clear();
                records.at(previous_vcf_record_position)->clear();
                merged_VCF.add_record_core(vcf_record_merged);
                previous_vcf_record_position = records.size() - 1;
                vcf_record_merged = *(records.at(previous_vcf_record_position));
                 */
            }
            vcf_record_merged.merge_record_into_this(vcf_record_to_be_merged_in);
        }else {
            merged_VCF.add_record_core(vcf_record_merged);
            vcf_record_merged = vcf_record_to_be_merged_in;
        }
    });
    merged_VCF.add_record_core(vcf_record_merged);

    merged_VCF.sort_records();
    assert(merged_VCF.get_VCF_size() <= vcf_size);

    return merged_VCF;
}

void VCF::correct_dot_alleles(const std::string &vcf_ref, const std::string &chrom) {
    //NB need to merge multiallelic before
    //NB cannot add covgs after
    auto vcf_size = records.size();
    for (auto &recordPointer : records) {
        auto &record = *recordPointer;
        if (record.chrom != chrom)
            continue;
        assert(vcf_ref.length() >= record.pos || assert_msg("vcf_ref.length() = " << vcf_ref.length() << "!>= record.pos "
                                                                                 << record.pos << "\n" << record.to_string(true, false) << "\n"
                                                                                 << vcf_ref));
        bool add_prev_letter = record.contains_dot_allele();

        if (add_prev_letter and record.pos > 0) {
            BOOST_LOG_TRIVIAL(debug) << record.pos;
            auto prev_letter = vcf_ref[record.pos - 1];
            BOOST_LOG_TRIVIAL(debug) << prev_letter;
            if (record.ref == "" or record.ref == ".")
                record.ref = prev_letter;
            else
                record.ref = prev_letter + record.ref;
                record.pos -= 1;
            for (auto &a : record.alts){
                if(a == "" or a == ".")
                    a = prev_letter;
                else
                    a = prev_letter + a;
            }


        } else if (add_prev_letter and record.pos + record.ref.length() + 1 < vcf_ref.length()) {
            auto next_letter = vcf_ref[record.pos + record.ref.length()];
            if (record.ref == "" or record.ref == ".") {
                next_letter = vcf_ref[record.pos];
                record.ref = next_letter;
            } else
                record.ref = record.ref + next_letter;
            for (auto &a : record.alts)
                if(a == "" or a == ".")
                    a = next_letter;
                else
                    a = a + next_letter;
        } else if (add_prev_letter) {
            record.clear();
        }
    }
    clean();
    assert(records.size() <= vcf_size);
    sort_records();
}

void VCF::make_gt_compatible() {
    //TODO: this can be a source of inneficiency, fix this?
    //TODO: interval tree does not to be indexed everytime this function is called, only if the interval tree changed
    BOOST_LOG_TRIVIAL(info) << now() << "Make all genotypes compatible";

    // set the GT compatible to the GT from coverage
    for (auto &recordPointer : records) {
        for (uint i = 0; i < recordPointer->sampleIndex_to_sampleInfo.size(); ++i) {
            auto gt_from_coverage = recordPointer->sampleIndex_to_sampleInfo[i].get_gt_coverages();
            recordPointer->sampleIndex_to_sampleInfo[i].set_gt_coverages_compatible(gt_from_coverage);
        }
    }


    // update GT compatible
    for (auto &pair : chrom2recordIntervalTree)
        pair.second.index();

    uint record_count=0;
    for (auto &recordPointer : records) {
        auto &record = *recordPointer;
        record_count++;
        for (uint i = 0; i < record.sampleIndex_to_sampleInfo.size(); ++i) {
            //retrieve all VCF records overlapping this record
            std::vector<VCFRecord*> overlappingVCFRecords = get_all_records_overlapping_the_given_record(record);

            for (VCFRecord* other_record_pointer : overlappingVCFRecords) {
                VCFRecord &other_record = *other_record_pointer;

                bool same_record_no_need_to_process = record == other_record;
                if (same_record_no_need_to_process)
                    continue;

                bool other_record_start_overlaps_this_record = record.pos <= other_record.pos && other_record.pos <= record.pos + record.ref.length();
                if (   other_record_start_overlaps_this_record
                    && record.sampleIndex_to_sampleInfo[i].is_gt_from_coverages_valid()
                    && other_record.sampleIndex_to_sampleInfo[i].is_gt_from_coverages_valid()) {

                    if (record.sampleIndex_to_sampleInfo[i].get_gt_coverages() == 0 and other_record.sampleIndex_to_sampleInfo[i].get_gt_coverages() == 0)
                        continue;
                    else if (record.sampleIndex_to_sampleInfo[i].get_likelihood_of_GT_from_coverages() >
                            other_record.sampleIndex_to_sampleInfo[i].get_likelihood_of_GT_from_coverages()) {
                        if (record.sampleIndex_to_sampleInfo[i].get_gt_coverages() == 0)
                            other_record.sampleIndex_to_sampleInfo[i].set_gt_coverages_compatible(0);
                        else
                            other_record.sampleIndex_to_sampleInfo[i].set_gt_coverages_compatible(boost::none);
                    } else {
                        if (other_record.sampleIndex_to_sampleInfo[i].get_gt_coverages() == 0)
                            record.sampleIndex_to_sampleInfo[i].set_gt_coverages_compatible(0);
                        else
                            record.sampleIndex_to_sampleInfo[i].set_gt_coverages_compatible(boost::none);
                    }
                }
            }
        }
    }
}

std::vector<VCFRecord*> VCF::get_all_records_overlapping_the_given_record (const VCFRecord &vcf_record) const {
    std::vector<VCFRecord*> overlappingVCFRecords;
    std::vector<size_t> overlaps;
    chrom2recordIntervalTree.at(vcf_record.chrom).overlap(vcf_record.pos, vcf_record.pos + vcf_record.ref.length() + 1, overlaps);
    for (size_t i = 0; i < overlaps.size(); ++i)
        overlappingVCFRecords.push_back(chrom2recordIntervalTree.at(vcf_record.chrom).data(overlaps[i]));

    //TODO: not sure if we need to sort, but doing it just in case... Check this later
    std::sort(overlappingVCFRecords.begin(), overlappingVCFRecords.end(), [](VCFRecord* lhs, VCFRecord* rhs) {
        return *lhs < *rhs;
    });

    return overlappingVCFRecords;
}



void VCF::save(const std::string &filepath, bool output_dot_allele, bool graph_is_simple, bool graph_is_nested, bool graph_has_too_many_alts, bool sv_type_is_snp, bool sv_type_is_indel,
               bool sv_type_is_ph_snps, bool sv_type_is_complex) {
    BOOST_LOG_TRIVIAL(debug) << "Saving VCF to " << filepath;
    std::ofstream handle;
    handle.open(filepath);
    handle << this->to_string(output_dot_allele, graph_is_simple, graph_is_nested, graph_has_too_many_alts, sv_type_is_snp, sv_type_is_indel,
            sv_type_is_ph_snps, sv_type_is_complex);
    handle.close();
    BOOST_LOG_TRIVIAL(debug) << "Finished saving " << this->records.size() << " entries to file";
}

std::string VCF::header() const {
    // find date
    time_t t = time(0);
    char mbstr[10];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));

    std::set<std::string> chroms;
    for (const auto record : records) {
        chroms.insert(record->chrom);
    }

    std::string header;
    header += "##fileformat=VCFv4.3\n";
    header += "##fileDate==";
    header += mbstr;
    header += "\n##ALT=<ID=SNP,Description=\"SNP\">\n";
    header += "##ALT=<ID=PH_SNPs,Description=\"Phased SNPs\">\n";
    header += "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n";
    header += "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and indels\">\n";
    header += "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n";
    header += "##ALT=<ID=SIMPLE,Description=\"Graph bubble is simple\">\n";
    header += "##ALT=<ID=NESTED,Description=\"Variation site was a nested feature in the graph\">\n";
    header += "##ALT=<ID=TOO_MANY_ALTS,Description=\"Variation site was a multinested feature with too many alts to include all in the VCF\">\n";
    header += "##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description=\"Type of graph feature\">\n";
    header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    header += "##FORMAT=<ID=MEAN_FWD_COVG,Number=A,Type=Integer,Description=\"Mean forward coverage\">\n";
    header += "##FORMAT=<ID=MEAN_REV_COVG,Number=A,Type=Integer,Description=\"Mean reverse coverage\">\n";
    header += "##FORMAT=<ID=MED_FWD_COVG,Number=A,Type=Integer,Description=\"Med forward coverage\">\n";
    header += "##FORMAT=<ID=MED_REV_COVG,Number=A,Type=Integer,Description=\"Med reverse coverage\">\n";
    header += "##FORMAT=<ID=SUM_FWD_COVG,Number=A,Type=Integer,Description=\"Sum forward coverage\">\n";
    header += "##FORMAT=<ID=SUM_REV_COVG,Number=A,Type=Integer,Description=\"Sum reverse coverage\">\n";
    header += "##FORMAT=<ID=GAPS,Number=A,Type=Float,Description=\"Number of gap bases\">\n";
    header += "##FORMAT=<ID=LIKELIHOOD,Number=A,Type=Float,Description=\"Likelihood\">\n";
    header += "##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description=\"Genotype confidence\">\n";
    for (const auto chrom : chroms){
        header += "##contig=<ID=" + chrom + ">\n";
    }
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (uint32_t i = 0; i != samples.size(); ++i) {
        header += "\t" + samples[i];
    }
    header += "\n";
    return header;
}

std::string VCF::to_string(bool genotyping_from_maximum_likelihood, bool genotyping_from_coverage,
                           bool output_dot_allele, bool graph_is_simple, bool graph_is_nested, bool graph_has_too_many_alts, bool sv_type_is_snp, bool sv_type_is_indel,
                           bool sv_type_is_ph_snps, bool sv_type_is_complex) {
    bool only_one_flag_is_set = ((int)(genotyping_from_maximum_likelihood) + (int)(genotyping_from_coverage)) == 1;
    assert(only_one_flag_is_set);


    std::stringstream out;
    out << header();

    // TODO: a side-effect of saving a VCF is sorting it, this might not be desirable
    // TODO: remove this side effect or always keep the VCF sorted
    sort_records();

    for (const auto &record : this->records) {
        bool record_has_dot_allele_and_should_be_output = output_dot_allele and record->contains_dot_allele();


        bool graph_type_condition_is_satisfied = (graph_is_simple and record->graph_type_is_simple()) or
                                                 (graph_is_nested and record->graph_type_is_nested()) or
                                                 (graph_has_too_many_alts and record->graph_type_has_too_many_alts());
        bool sv_type_condition_is_satisfied = (sv_type_is_snp and record->svtype_is_SNP()) or
                                              (sv_type_is_indel and record->svtype_is_indel()) or
                                              (sv_type_is_ph_snps and record->svtype_is_PH_SNPs()) or
                                              (sv_type_is_complex and record->svtype_is_complex());
        bool graph_and_sv_type_conditions_are_satisfied = graph_type_condition_is_satisfied and sv_type_condition_is_satisfied;

        bool record_should_be_output = record_has_dot_allele_and_should_be_output or graph_and_sv_type_conditions_are_satisfied;

        if (record_should_be_output) {
            out << record->to_string(genotyping_from_maximum_likelihood, genotyping_from_coverage) << std::endl;
        }
    }

    return out.str();
}


// TODO: check if we keep this, it is only used in tests - better to keep in a VCFMock class
/*
void VCF::load(const std::string &filepath) {
    BOOST_LOG_TRIVIAL(debug) << "Loading VCF from " << filepath;
    VCFRecord vr;
    std::string line;
    std::stringstream ss;
    uint32_t added = 0;
    std::vector<std::string> sample_names = {};
    // NB this doesn't currently clear records first. Do we want to?

    std::ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line[0] != '#') {
                ss << line;
                ss >> vr;
                ss.clear();
                add_record(vr, sample_names);
                added += 1;
            } else if (line[1] != '#'){
                auto sample_string = line.replace(0,45, "");
                sample_names = split(sample_string, "\t");
            }
        }
    } else {
        std::cerr << "Unable to open VCF file " << filepath << std::endl;
        std::exit(1);
    }
    BOOST_LOG_TRIVIAL(debug) << "Finished loading " << added << " entries to VCF, which now has size "
                             << records.size();
}
 */

bool VCF::operator==(const VCF &y) const {
    if (records.size() != y.records.size()) { return false; }
    for (uint32_t i = 0; i != y.records.size(); ++i) {
        if (find_record_in_records(*(y.records[i])) == records.end()) {
            return false;
        }
    }
    return true;
}

bool VCF::operator!=(const VCF &y) const {
    return !(*this == y);
}


void VCF::concatenateVCFs(const std::vector<std::string> &VCFPathsToBeConcatenated, const std::string &sink) {
    std::ofstream outFile(sink);
    bool headerIsOutput = false;
    for (const auto &VCFPath : VCFPathsToBeConcatenated) {
        std::ifstream inFile(VCFPath);
        std::string line;
        while (std::getline(inFile, line)) {
            if ((line[0]!='#') || (line[0]=='#' && headerIsOutput==false))
                outFile << line << std::endl;
        }
        headerIsOutput = true;
        inFile.close();
    }
    outFile.close();
}
