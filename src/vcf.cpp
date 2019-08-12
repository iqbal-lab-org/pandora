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
    std::unordered_map<std::string, std::vector<uint16_t>> empty_map;
    VCFRecord vr(c, p, r, a, i, g);
    if (find_record_in_records(vr) == records.end()) { //TODO: improve this search to log(n) using a map or sth
        add_record_core(vr);
        records.back()->samples.insert(records.back()->samples.end(), samples.size(), empty_map);
    }
    //cout << "added record: " << vr << endl;
}



VCFRecord &VCF::add_record(VCFRecord &vr, const std::vector<std::string> &sample_names) {
    assert(vr.samples.size() == sample_names.size() or sample_names.size() == 0);

    std::unordered_map<std::string, std::vector<uint16_t>> empty_map;
    auto record_it = find_record_in_records(vr); //TODO: improve this search to log(n) using a map or sth
    if (record_it == records.end()) {
        add_record_core(vr);
        records.back()->samples.clear();
        records.back()->samples.insert(records.back()->samples.end(), samples.size(), empty_map);
        record_it = --records.end();
    }

    for (uint32_t i=0; i < sample_names.size(); ++i){
        auto &name = sample_names[i];
        ptrdiff_t sample_index = get_sample_index(name);
        (*record_it)->samples[sample_index] = vr.samples[i];
    }

    return **record_it;
}

void VCF::add_samples(const std::vector<std::string> sample_names) {
    for (uint32_t i=0; i < sample_names.size(); ++i){
        auto &name = sample_names[i];
        ptrdiff_t sample_index = get_sample_index(name);
    }
}

void VCF::add_formats(const std::vector<std::string> &v) {
    for (auto &record : records) {
        record->add_formats(v);
    }
}

ptrdiff_t VCF::get_sample_index(const std::string &name) {
    std::unordered_map<std::string, std::vector<uint16_t>> empty_map;

    // if this sample has not been added before, add a column for it
    auto sample_it = find(samples.begin(), samples.end(), name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(name);
        for (uint32_t i = 0; i != records.size(); ++i) {
            records[i]->samples.push_back(empty_map);
            assert(samples.size() == records[i]->samples.size());
        }
        return samples.size() - 1;
    } else {
        return distance(samples.begin(), sample_it);
    }
}

void VCF::add_sample_gt(const std::string &name, const std::string &c, const uint32_t p, const std::string &r,
                        const std::string &a) {
    //cout << "adding gt " << c << " " << p << " " << r << " vs " << name << " " << a << endl;
    if (r == "" and a == "") {
        return;
    }

    //cout << "adding gt " << r << " vs " << a << endl;

    ptrdiff_t sample_index = get_sample_index(name);

    VCFRecord vr(c, p, r, a);
    VCFRecord *vrp;
    bool added = false;
    auto it = find_record_in_records(vr); //TODO: improve this search to log(n) using a map or sth
    if (it != records.end()) {
        //cout << "found record with this ref and alt" << endl;
        (*it)->samples[sample_index]["GT"] = {1};
        vrp = it->get();
    } else {
        //cout << "didn't find a record for pos " << p << " ref " << r << " and alt " << a << endl;
        // either we have the ref allele, an alternative allele for a too nested site, or a mistake
        for (uint32_t i = 0; i != records.size(); ++i) {
            if (records[i]->chrom == c and records[i]->pos == p and r == a and records[i]->ref == r) {
                //cout << "have ref allele" << endl;
                records[i]->samples[sample_index]["GT"] = {0};
                vrp = records[i].get();
                added = true;
            }
        }
        if (!added and r != a) {
            //cout << "have new allele not in graph" << endl;
            add_record(c, p, r, a, "SVTYPE=COMPLEX", "GRAPHTYPE=TOO_MANY_ALTS");
            records.back()->samples[sample_index]["GT"] = {1};
            vrp = records.back().get();
            added = true;
        }
        // check not mistake
        assert(added);
    }

    // update other samples at this site if they have ref allele at this pos
    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i]->chrom == c and records[i]->pos <= p and records[i]->pos + records[i]->ref.length() > p) {
            for (uint32_t j = 0; j != records[i]->samples.size(); ++j) {
                if (records[i]->samples[j].find("GT") != records[i]->samples[j].end()
                    and !records[i]->samples[j]["GT"].empty()
                    and records[i]->samples[j]["GT"][0] == 0) {
                    //cout << "update my record to have ref allele also for sample " << j << endl;
                    //cout << records[i] << endl;
                    vrp->samples[j]["GT"] = {0};
                }
            }
        }
    }
}

void VCF::add_sample_ref_alleles(const std::string &sample_name, const std::string &chrom, const uint32_t &pos,
                                 const uint32_t &pos_to) {

    ptrdiff_t sample_index = get_sample_index(sample_name);

    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i]->chrom == chrom and pos <= records[i]->pos and
            records[i]->pos + records[i]->ref.length() <= pos_to) {
            records[i]->samples[sample_index]["GT"] = {0};
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
    std::unordered_map<std::string, std::vector<uint16_t>> empty_u_map;
    assert(original_size < std::numeric_limits<uint_least64_t>::max() ||
           assert_msg("VCF size has got too big to use the append feature"));
    for (uint_least64_t i = 0; i < original_size; ++i) {
        records[i]->samples.insert(records[i]->samples.end(), num_samples_added, empty_u_map);
    }

    BOOST_LOG_TRIVIAL(debug) << "add the " << other_vcf.records.size() << " records";
    for (auto &recordPointer : other_vcf.records) {
        auto &record = *recordPointer;
        VCFRecord &vr = add_record(record, other_vcf.samples);
        for (uint_least16_t j = 0; j < other_vcf.samples.size(); ++j) {
            vr.samples[other_sample_positions[j]] = record.samples[j];
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

void VCF::genotype(const std::vector<uint32_t> &expected_depth_covg, const float &error_rate,
                   const uint16_t confidence_threshold,
                   const uint32_t &min_allele_covg, const float &min_fraction_allele_covg,
                   const uint32_t &min_site_total_covg, const uint32_t &min_site_diff_covg, bool snps_only) {
    BOOST_LOG_TRIVIAL(info) << now() << "Genotype VCF";
    for (auto &vrPointer : records) {
        auto &vr = *vrPointer;
        if (not snps_only or (vr.ref.length() == 1 and !vr.alt.empty() and vr.alt[0].length() == 1)) {
            vr.likelihood(expected_depth_covg, error_rate, min_allele_covg, min_fraction_allele_covg);
            vr.confidence(min_site_total_covg, min_site_diff_covg);
            vr.genotype(confidence_threshold);
        }
    }
    add_formats({"GT_CONF", "LIKELIHOOD"});
    BOOST_LOG_TRIVIAL(info) << now() << "Make all genotypes compatible";
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

void merge_sample_key(std::unordered_map<std::string, std::vector<uint16_t>> &first,
                      const std::unordered_map<std::string, std::vector<uint16_t>> &second,
                      const std::string &key) {
    if (first.empty() or second.empty() or first.find(key) == first.end() or first[key].empty()) {
        return;
    } else if (first.find(key) != first.end() and (second.find(key) == second.end() or second.at(key).empty())) {
        first.erase(key);
    } else if (first[key][0] == second.at(key).at(0)) {
        bool ref = true;
        for (const auto &val : second.at(key)) {
            if (!ref)
                first[key].push_back(val);
            ref = false;
        }
    } else
        first.erase(key);
}

void merge_regt_sample_key(std::unordered_map<std::string, std::vector<float>> &first,
                           const std::unordered_map<std::string, std::vector<float>> &second,
                           const std::string &key) {
    if (first.empty() or second.empty() or first.find(key) == first.end() or first[key].empty()) {
        return;
    } else if (first.find(key) != first.end() and (second.find(key) == second.end() or second.at(key).empty())) {
        first.erase(key);
    } else if (first[key][0] == second.at(key).at(0)) {
        bool ref = true;
        for (const auto &val : second.at(key)) {
            if (!ref)
                first[key].push_back(val);
            ref = false;
        }
    } else
        first.erase(key);
}

void merge_gt(VCFRecord &first, const VCFRecord &second, const uint16_t i, const uint16_t prev_alt_size) {
    if (first.samples.size() < i or second.samples.size() < i) {
        return;
    } else if (second.samples[i].find("GT") == second.samples[i].end()
               or second.samples[i].at("GT").empty()) {
        return;
    } else if (first.samples[i].find("GT") == first.samples[i].end()
               or first.samples[i]["GT"].empty()) {
        if (second.samples[i].at("GT")[0] == 0) {
            first.samples[i]["GT"] = {0};
        } else {
            uint16_t new_allele = second.samples[i].at("GT")[0] + prev_alt_size;
            first.samples[i]["GT"] = {new_allele};
        }
    } else if (first.samples[i]["GT"][0] != 0 or second.samples[i].at("GT")[0] != 0) {
        //conflict, try to resolve with likelihoods
        if (first.regt_samples.size() > i
            and first.regt_samples[i].find("LIKELIHOOD") != first.regt_samples[i].end()) {
            first.confidence();
            first.genotype(5);
        } else {
            first.samples[i]["GT"] = {};
        }
    }
}


void VCF::merge_multi_allelic(uint32_t max_allele_length) {
    if (records.size() < 2)
        return;

    uint32_t prev_pos = 0;
    VCFRecord prev_vr(*(records.at(prev_pos)));
    auto vcf_size = records.size();
    auto reserve_size = vcf_size * 1.05;
    records.reserve(reserve_size);
    for (uint32_t current_pos = 1; current_pos < vcf_size; ++current_pos) {
        const auto &record = *(records.at(current_pos));
        //cout << "comparing record " << current_pos << "/" << vcf_size << " to record " << prev_pos << endl;

        if (record != prev_vr
            and prev_vr.chrom == record.chrom
            and prev_vr.pos == record.pos
            and prev_vr.ref == record.ref
            and prev_vr.ref != "."
            and prev_vr.ref != ""
            and prev_vr.ref.length() <= max_allele_length
            and prev_vr.alt[0].length() <= max_allele_length) {

            // merge alts
            uint16_t prev_alt_size = prev_vr.alt.size();
            bool short_enough = true;
            for (const auto &a : record.alt) {
                if (a.length() > max_allele_length)
                    short_enough = false;
                prev_vr.alt.push_back(a);
            }
            if (!short_enough) {
                prev_pos = current_pos;
                prev_vr = *(records.at(prev_pos));
                continue;
            }

            // merge count/likelihood data
            if (record.samples.empty()) {
                records.at(current_pos)->clear();
                records.at(prev_pos)->clear();
                add_record_core(prev_vr);
                prev_pos = records.size() - 1;
                prev_vr = *(records.at(prev_pos));
            }
            for (uint i = 0; i < record.samples.size(); ++i) {
                auto keys = {"MEAN_FWD_COVG", "MEAN_REV_COVG",
                             "MED_FWD_COVG", "MED_REV_COVG",
                             "SUM_FWD_COVG", "SUM_REV_COVG"};
                for (const auto &key: keys) {
                    merge_sample_key(prev_vr.samples[i], record.samples[i], key);
                }
                keys = {"LIKELIHOOD", "GT_CONF", "GAPS"};
                if (!prev_vr.regt_samples.empty() and !record.regt_samples.empty()) {
                    for (const auto &key: keys)
                        merge_regt_sample_key(prev_vr.regt_samples[i], record.regt_samples[i], key);
                }

                //merge GT
                merge_gt(prev_vr, record, i, prev_alt_size);
                records.at(current_pos)->clear_sample(i);
                records.at(prev_pos)->clear_sample(i);
            }
            add_record_core(prev_vr);
            prev_pos = records.size() - 1;
            prev_vr = *(records.at(prev_pos));
        } else if (record != prev_vr) {
            prev_pos = current_pos;
            prev_vr = *(records.at(prev_pos));
        }
    }
    clean();
    assert(records.size() <= vcf_size);
    sort_records();
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
                                                                                 << record.pos << "\n" << record << "\n"
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
            for (auto &a : record.alt){
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
            for (auto &a : record.alt)
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
    for (auto &pair : chrom2recordIntervalTree)
        pair.second.index();

    uint record_count=0;
    for (auto &recordPointer : records) { //goes through all records
        auto &record = *recordPointer;
        record_count++;
        for (uint i = 0; i < record.samples.size(); ++i) { //goes through all samples of this record

            //retrieve all VCF records overlapping this record
            std::vector<VCFRecord*> overlappingVCFRecords;

            std::vector<size_t> overlaps;
            chrom2recordIntervalTree[record.chrom].overlap(record.pos, record.pos + record.ref.length() + 1, overlaps);
            for (size_t i = 0; i < overlaps.size(); ++i)
                overlappingVCFRecords.push_back(chrom2recordIntervalTree[record.chrom].data(overlaps[i]));

            //TODO: not sure if we need to sort, but doing it just in case... Check this later
            std::sort(overlappingVCFRecords.begin(), overlappingVCFRecords.end(), [](VCFRecord* lhs, VCFRecord* rhs) {
                return *lhs < *rhs;
            });

            for (VCFRecord* other_record_pointer : overlappingVCFRecords) {
                VCFRecord &other_record = *other_record_pointer;
                if (record == other_record) //no need to process this
                    continue;
                if (   record.pos <= other_record.pos
                    && other_record.pos <= record.pos + record.ref.length()
                    && record.samples[i].find("GT") != record.samples[i].end()
                    && other_record.samples[i].find("GT") != other_record.samples[i].end()
                    && record.samples[i]["GT"].size() > 0
                    && other_record.samples[i]["GT"].size() > 0) {

                    if (record.samples[i]["GT"][0] == 0 and other_record.samples[i]["GT"][0] == 0)
                        continue;
                    else if (not record.regt_samples.empty() and not other_record.regt_samples.empty()
                             and record.regt_samples[i].find("LIKELIHOOD") != record.regt_samples[i].end()
                             and
                             other_record.regt_samples[i].find("LIKELIHOOD") != other_record.regt_samples[i].end()) {
                        if (record.regt_samples[i]["LIKELIHOOD"][record.samples[i]["GT"][0]] >
                            other_record.regt_samples[i]["LIKELIHOOD"][other_record.samples[i]["GT"][0]]) {
                            if (record.samples[i]["GT"][0] == 0)
                                other_record.samples[i]["GT"] = {0};
                            else
                                other_record.samples[i]["GT"] = {};
                        } else {
                            if (other_record.samples[i]["GT"][0] == 0)
                                record.samples[i]["GT"] = {0};
                            else
                                record.samples[i]["GT"] = {};
                        }
                    } else {
                        other_record.samples[i] = {};
                        record.samples[i] = {};
                    }
                }
            }
        }
    }
}

std::string VCF::header() {
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
    std::cout << header << std::endl;
    return header;
}

// NB in the absence of filter flags being set to true, all results are saved. If one or more filter flags for SVTYPE are set, 
// then only those matching the filter are saved. Similarly for GRAPHTYPE.
void VCF::save(const std::string &filepath, bool simple, bool complexgraph, bool toomanyalts, bool snp, bool indel,
               bool phsnps, bool complexvar) {
    /*if (samples.empty())
    {
	    cout << now() << "Did not save VCF for sample" << endl;
	    return;
    }*/
    BOOST_LOG_TRIVIAL(debug) << "Saving VCF to " << filepath;

    // open and write header
    std::ofstream handle;
    handle.open(filepath);

    handle << header();

    sort_records();

    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i]->contains_dot_allele())
            continue;

        if (((!simple and !complexgraph) or
             (simple and records[i]->info.find("GRAPHTYPE=SIMPLE") != std::string::npos) or
             (complexgraph and records[i]->info.find("GRAPHTYPE=NESTED") != std::string::npos) or
             (toomanyalts and records[i]->info.find("GRAPHTYPE=TOO_MANY_ALTS") != std::string::npos)) and
            ((!snp and !indel and !phsnps and !complexvar) or
             (snp and records[i]->info.find("SVTYPE=SNP") != std::string::npos) or
             (indel and records[i]->info.find("SVTYPE=INDEL") != std::string::npos) or
             (phsnps and records[i]->info.find("SVTYPE=PH_SNPs") != std::string::npos) or
             (complexvar and records[i]->info.find("SVTYPE=COMPLEX") != std::string::npos))) {
            handle << *(records[i]);
        }
    }
    handle.close();
    BOOST_LOG_TRIVIAL(debug) << "Finished saving " << records.size() << " entries to file";
}

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

std::ostream &operator<<(std::ostream &out, VCF const &m) {
    for (const auto &record : m.records) {
        out << (*record);
    }
    return out;
}
