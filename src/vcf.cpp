#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <vector>
#include <algorithm>
#include "vcfrecord.h"
#include "vcf.h"
#include "utils.h"
#include "localnode.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

VCF::VCF() {};

VCF::~VCF() {
    clear();
};

void VCF::add_record(string c, uint32_t p, string r, string a, string i, string g) {
    unordered_map<string,uint8_t> empty_map;
    VCFRecord vr(c, p, r, a, i, g);
    if (find(records.begin(), records.end(), vr) == records.end()) {
        records.push_back(vr);
        records.back().samples.insert(records.back().samples.end(), samples.size(), empty_map);
    }
    //cout << "added record: " << vr << endl;
}

VCFRecord& VCF::add_record(VCFRecord &vr) {
    unordered_map<string,uint8_t> empty_map;
    auto record_it = find(records.begin(), records.end(), vr);
    if (record_it == records.end()) {
        records.push_back(vr);
        if (samples.size() > vr.samples.size())
            records.back().samples.insert(records.back().samples.end(), samples.size()-vr.samples.size(), empty_map);
        return records.back();
    }
    return *record_it;
    // handle mismatched sample data?
}

void VCF::add_formats(const vector<string>& v){
    for(auto record : records){
        record.add_formats(v);
    }
}

ptrdiff_t VCF::get_sample_index(const string& name){
    unordered_map<string,uint8_t> empty_map;

    // if this sample has not been added before, add a column for it
    vector<string>::iterator sample_it = find(samples.begin(), samples.end(), name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(name);
        for (uint32_t i = 0; i != records.size(); ++i) {
            records[i].samples.push_back(empty_map);
            assert(samples.size() == records[i].samples.size());
        }
        return samples.size() - 1;
    } else {
        return distance(samples.begin(), sample_it);
    }
}

void VCF::add_sample_gt(const string &name, const string &c, const uint32_t p, const string &r, const string &a) {
    //cout << "adding gt " << c << " " << p << " " << r << " vs " << name << " " << a << endl;
    if (r == "" and a == "") {
        return;
    }

    //cout << "adding gt " << r << " vs " << a << endl;

    ptrdiff_t sample_index = get_sample_index(name);

    VCFRecord vr(c, p, r, a);
    VCFRecord* vrp;
    bool added = false;
    vector<VCFRecord>::iterator it = find(records.begin(), records.end(), vr);
    if (it != records.end()) {
        //cout << "found record with this ref and alt" << endl;
        it->samples[sample_index]["GT"] = 1;
        vrp = &*it;
    } else {
        //cout << "didn't find a record for pos " << p << " ref " << r << " and alt " << a << endl;
        // either we have the ref allele, an alternative allele for a too nested site, or a mistake
        for (uint32_t i = 0; i != records.size(); ++i) {
            if (records[i].chrom == c and records[i].pos == p and r == a and records[i].ref == r) {
                //cout << "have ref allele" << endl;
                records[i].samples[sample_index]["GT"] = 0;
                vrp = &records[i];
                added = true;
            }
        }
        if (added == false and r != a) {
            //cout << "have new allele not in graph" << endl;
            add_record(c, p, r, a, "SVTYPE=COMPLEX", "GRAPHTYPE=TOO_MANY_ALTS");
            records.back().samples[sample_index]["GT"] = 1;
            vrp = &records.back();
            added = true;
        }
        // check not mistake
        assert(added == true);
    }

    // update other samples at this site if they have ref allele at this pos
    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i].chrom == c and records[i].pos <= p and records[i].pos + records[i].ref.length() > p) {
            for (uint32_t j = 0; j != records[i].samples.size(); ++j) {
                if (records[i].samples[j].find("GT") != records[i].samples[j].end() and records[i].samples[j]["GT"] == 0) {
                    //cout << "update my record to have ref allele also for sample " << j << endl;
                    //cout << records[i] << endl;
                    vrp->samples[j]["GT"] = 0;
                }
            }
        }
    }

    return;
}

void VCF::add_sample_ref_alleles(const string &sample_name, const string &chrom, const uint32_t &pos, const uint32_t &pos_to) {
    ptrdiff_t sample_index;
    unordered_map<string,uint8_t> empty_map;

    // if this sample has not been added before, add a column for it
    vector<string>::iterator sample_it = find(samples.begin(), samples.end(), sample_name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(sample_name);
        for (uint32_t i = 0; i != records.size(); ++i) {
            records[i].samples.push_back(empty_map);
        }
        sample_index = samples.size() - 1;
    } else {
        sample_index = distance(samples.begin(), sample_it);
    }

    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i].chrom == chrom and pos <= records[i].pos and records[i].pos + records[i].ref.length() <= pos_to) {
            records[i].samples[sample_index]["GT"] = 0;
            //cout << "update record " << records[i] << endl;
        }
    }
    return;
}

void VCF::clear() {
    records.clear();
}

void VCF::append_vcf(const VCF &other_vcf){
    auto original_size = records.size();
    auto num_samples_added = 0;

    vector<uint_least16_t > other_sample_positions;
    for (const auto sample : other_vcf.samples){
        auto sample_it = find(samples.begin(), samples.end(), sample);
        if (sample_it == samples.end()){
            samples.push_back(sample);
            other_sample_positions.push_back(samples.size() - 1);
            num_samples_added += 1;
        } else {
            other_sample_positions.push_back(distance(samples.begin(), sample_it));
        }
    }

    unordered_map<string,uint8_t> empty_u_map;
    for (uint_least16_t i=0; i<original_size; ++i){
        records[i].samples.insert(records[i].samples.end(), num_samples_added, empty_u_map);
    }

    for (auto record : other_vcf.records){
        VCFRecord& vr = add_record(record);
        for(uint_least16_t j=0; j<other_vcf.samples.size(); ++j){
            vr.samples[other_sample_positions[j]] = record.samples[j];
            //NB this overwrites old data without checking
        }
    }
}

void VCF::sort_records() {
    sort(records.begin(), records.end());
    return;
}

bool VCF::pos_in_range(const uint32_t from, const uint32_t to, const string& chrom) const {
    // is there a record contained in the range from,to?
    for (auto record : records) {
        if (chrom == record.chrom and from < record.pos and record.pos + record.ref.length() <= to) {
            return true;
        }
    }
    return false;
}

void VCF::regenotype(const uint32_t & expected_depth_covg, const float & error_rate, const uint8_t confidence_threshold) {
    for (auto &vr : records){
        if (vr.ref.length() == 1 and vr.alt.length() == 1)
        {
            vr.likelihood(expected_depth_covg,error_rate);
            vr.confidence();
            vr.regenotype(confidence_threshold);
        }
    }
    add_formats({"GT_CONF", "REF_LIKELIHOOD", "ALT_LIKELIHOOD"});
}

string VCF::header() {
    // find date
    time_t t = time(0);
    char mbstr[10];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));

    string header;
    header += "##fileformat=VCFv4.3\n";
    header += "##fileDate==";
    header += mbstr;
    header += "\n##ALT=<ID=SNP,Description=\"SNP\">\n" ;
    header += "##ALT=<ID=PH_SNPs,Description=\"Phased SNPs\">\n";
    header += "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n";
    header += "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and indels\">\n";
    header += "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n";
    header += "##ALT=<ID=SIMPLE,Description=\"Graph bubble is simple\">\n";
    header += "##ALT=<ID=NESTED,Description=\"Variation site was a nested feature in the graph\">\n";
    header += "##ALT=<ID=TOO_MANY_ALTS,Description=\"Variation site was a multinested feature with too many alts to include all in the VCF\">\n";
    header += "##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description=\"Type of graph feature\">\n";
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (uint32_t i = 0; i != samples.size(); ++i) {
        header += "\t" + samples[i];
    }
    header += "\n";
    return header;
}

// NB in the absence of filter flags being set to true, all results are saved. If one or more filter flags for SVTYPE are set, 
// then only those matching the filter are saved. Similarly for GRAPHTYPE.
void VCF::save(const string &filepath, bool simple, bool complexgraph, bool toomanyalts, bool snp, bool indel,
               bool phsnps, bool complexvar) {
    /*if (samples.size() == 0)
    {
	    cout << now() << "Did not save VCF for sample" << endl;
	    return;
    }*/
    cout << now() << "Saving VCF to " << filepath << endl;

    // open and write header
    ofstream handle;
    handle.open(filepath);

    handle << header();

    sort_records();

    for (uint32_t i = 0; i != records.size(); ++i) {
        if (((simple == false and complexgraph == false) or
             (simple == true and records[i].info.find("GRAPHTYPE=SIMPLE") != std::string::npos) or
             (complexgraph == true and records[i].info.find("GRAPHTYPE=NESTED") != std::string::npos) or
             (toomanyalts == true and records[i].info.find("GRAPHTYPE=TOO_MANY_ALTS") != std::string::npos)) and
            ((snp == false and indel == false and phsnps == false and complexvar == false) or
             (snp == true and records[i].info.find("SVTYPE=SNP") != std::string::npos) or
             (indel == true and records[i].info.find("SVTYPE=INDEL") != std::string::npos) or
             (phsnps == true and records[i].info.find("SVTYPE=PH_SNPs") != std::string::npos) or
             (complexvar == true and records[i].info.find("SVTYPE=COMPLEX") != std::string::npos))) {
            handle << records[i];
        }
    }
    handle.close();
    cout << now() << "Finished saving " << records.size() << " entries to file" << endl;
    return;
}

void VCF::load(const string &filepath) {
    cout << now() << "Loading VCF from " << filepath << endl;
    VCFRecord vr;
    string line;
    stringstream ss;
    uint32_t added = 0;
    // NB this doesn't currently clear records first. Do we want to?

    ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line[0] != '#') {
                ss << line;
                ss >> vr;
                ss.clear();
                add_record(vr);
                added += 1;
            }
        }
    } else {
        cerr << "Unable to open VCF file " << filepath << endl;
        exit(1);
    }
    cout << now() << "Finished loading " << added << " entries to VCF, which now has size " << records.size() << endl;
    return;
}

void VCF::write_aligned_fasta(const string &filepath, const string &chrom, const vector<LocalNodePtr> &lmp) {
    cout << now() << "Write aligned fasta to " << filepath << endl;
    sort_records();

    if (lmp.empty() or samples.empty()) {
        return;
    }

    /*cout << "lmp:" << endl;
    for (uint i=0; i!= lmp.size(); ++i)
    {
	cout << *lmp[i] << endl;
    }*/

    int prev_pos = -1;
    uint32_t max_len = 0;
    uint32_t n = 0; // position in lmp
    uint32_t ref_len = 0;
    vector<string> seqs(samples.size(), "");
    vector<uint32_t> alt_until(samples.size(), 0);

    for (uint32_t i = 0; i != records.size(); ++i) {
        if (records[i].chrom == chrom and records[i].pos != (uint32_t) prev_pos) {
            // equalise lengths of sequences
            for (uint32_t j = 0; j != samples.size(); ++j) {
                if (seqs[j].length() < max_len) {
                    //cout << "equalise length of seq " << j << endl;
                    string s(max_len - seqs[j].length(), '-'); // s == "------"
                    seqs[j] += s;
                }
            }

            // add ref sequence for gaps
            while (records[i].chrom == chrom and ref_len < records[i].pos and n < lmp.size()) {
                for (uint32_t j = 0; j != samples.size(); ++j) {
                    if (alt_until[j] < records[i].pos) {
                        //cout << "add ref gap allele to seq " << j << endl;
                        seqs[j] += lmp[n]->seq;
                    }
                }
                ref_len += lmp[n]->seq.length();
                n++;
            }
        }
        //cout << "record[" << i << "] " << records[i] << endl;
        for (uint32_t j = 0; j != samples.size(); ++j) {
            if (records[i].chrom == chrom and records[i].samples[j]["GT"] == 0 and records[i].pos != (uint32_t) prev_pos and
                pos_in_range(records[i].pos, records[i].pos + records[i].ref.length(),chrom) == false) {
                //cout << j << " add ref allele at site" << endl;
                seqs[j] += records[i].ref;
                max_len = max(max_len, (uint32_t) seqs[j].length());
                ref_len += records[i].ref.length();
                n++;
            } else if (records[i].chrom == chrom and records[i].samples[j]["GT"] == 1) {
                //cout << j << " add alt allele at site" << endl;
                assert(records[i].pos != (uint32_t) prev_pos);
                seqs[j] += records[i].alt;
                max_len = max(max_len, (uint32_t) seqs[j].length());
                alt_until[j] = records[i].pos + records[i].ref.length();
            }
        }
    }

    // equalise lengths of sequences
    for (uint32_t j = 0; j != samples.size(); ++j) {
        if (seqs[j].length() < max_len) {
            //cout << "equalise length of seq " << j << endl;
            string s(max_len - seqs[j].length(), '-'); // s == "------"
            seqs[j] += s;
        }
    }

    // add ref sequence for end gaps
    while (n < lmp.size()) {
        for (uint32_t j = 0; j != samples.size(); ++j) {
            if (alt_until[j] <= ref_len) {
                //cout << "add ref gap allele to seq " << j << endl;
                seqs[j] += lmp[n]->seq;
                //} else {
                //	cout << "alt_until[j] > ref_len" << alt_until[j]  << " " << ref_len << endl;
            }
        }
        ref_len += lmp[n]->seq.length();
        n++;
    }
    /*cout << "seqs now are " << endl;
    for (uint j=0; j!=samples.size(); ++j)
    {
	cout << j << " " << seqs[j] << endl;
    }*/

    ofstream handle;
    handle.open(filepath);
    assert (!handle.fail());

    for (uint32_t j = 0; j != samples.size(); ++j) {
        handle << ">" << samples[j] << endl;
        handle << seqs[j] << endl;
    }
    handle.close();
    return;

}

bool VCF::operator==(const VCF &y) const {
    if (records.size() != y.records.size()) { return false; }
    for (uint32_t i = 0; i != y.records.size(); ++i) {
        if (find(records.begin(), records.end(), y.records[i]) == records.end()) {
            return false;
        }
    }
    return true;
}

bool VCF::operator!=(const VCF &y) const {
    return !(*this==y);
}

std::ostream &operator<<(std::ostream &out, VCF const &m) {
    for (auto record : m.records){
        out << record;
    }
    return out;
}
