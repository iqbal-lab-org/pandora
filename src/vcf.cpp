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
    VCFRecord vr(c, p, r, a, i, g);
    if (find(records.begin(), records.end(), vr) == records.end()) {
        records.push_back(vr);
        records.back().samples.insert(records.back().samples.end(), samples.size(), ".");
    }
    //cout << "added record: " << vr << endl;
}

void VCF::add_record(VCFRecord &vr) {
    if (find(records.begin(), records.end(), vr) == records.end()) {
        records.push_back(vr);
        records.back().samples.insert(records.back().samples.end(), samples.size(), ".");
    }
}

void VCF::add_sample_gt(const string &name, const string &c, const uint32_t p, const string &r, const string &a) {
    //cout << "adding gt " << c << " " << p << " " << r << " vs " << name << " " << a << endl;
    if (r == "" and a == "") {
        return;
    }

    //cout << "adding gt " << r << " vs " << a << endl;

    ptrdiff_t sample_index;

    // if this sample has not been added before, add a column for it
    vector<string>::iterator sample_it = find(samples.begin(), samples.end(), name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(name);
        for (uint i = 0; i != records.size(); ++i) {
            records[i].samples.push_back(".");
        }
        sample_index = samples.size() - 1;
    } else {
        sample_index = distance(samples.begin(), sample_it);
    }

    VCFRecord vr(c, p, r, a);
    bool added = false;
    vector<VCFRecord>::iterator it = find(records.begin(), records.end(), vr);
    if (it != records.end()) {
        it->samples[sample_index] = "1";
        //cout << "found record with this ref and alt" << endl;
        for (uint i = 0; i != records.size(); ++i) {
            if (records[i].pos == p and records[i].alt != a) {
                //assert(records[i].ref == r or r == "");
                records[i].samples[sample_index] = ".";
            } else if (records[i].pos > p) {
                break;
            }
        }
    } else {
        //cout << "didn't find a record for pos " << p << " ref " << r << " and alt " << a << endl;
        // either we have the ref allele, an alternative allele for a too nested site, or a mistake
        for (uint i = 0; i != records.size(); ++i) {
            if (records[i].pos == p and r == a and records[i].ref == r) {
                //cout << "have ref allele" << endl;
                //assert(records[i].ref == r or r == "" ||
                //assert_msg("at pos " << records[i].pos << " existing ref is " << records[i].ref << " which is not equal to " << r));
                records[i].samples[sample_index] = "0";
                added = true;
            } else if (records[i].pos == p and r != a) {
                //cout << "found another alt at the position" << endl;
                //assert(records[i].ref == r or r == "" || assert_msg("at pos " << records[i].pos << " existing ref is " << records[i].ref << " which is not equal to " << r));
                records[i].samples[sample_index] = ".";
                /*} else if (records[i].pos > p) {
		    break;*/
            }
        }
        if (added == false and r != a) {
            //cout << "have very nested allele" << endl;
            add_record(c, p, r, a, "SVTYPE=COMPLEX", "GRAPHTYPE=TOO_MANY_ALTS");
            records.back().samples[sample_index] = "1";
            added = true;
            // also check if other samples had ref allele at this pos
            for (uint i = 0; i != records.size(); ++i) {
                if (records[i].pos <= p and records[i].pos + records[i].ref.length() > p) {
                    for (uint j = 0; j != records[i].samples.size(); ++j) {
                        if (records[i].samples[j] == "0") {
                            records.back().samples[j] = "0";
                        }
                    }
                }
            }
        }
        assert(added == true);
    }

    return;
}

void VCF::add_sample_ref_alleles(const string &sample_name, const string &name, const uint &pos, const uint &pos_to) {
    ptrdiff_t sample_index;

    // if this sample has not been added before, add a column for it
    vector<string>::iterator sample_it = find(samples.begin(), samples.end(), sample_name);
    if (sample_it == samples.end()) {
        //cout << "this is the first time this sample has been added" << endl;
        samples.push_back(sample_name);
        for (uint i = 0; i != records.size(); ++i) {
            records[i].samples.push_back(".");
        }
        sample_index = samples.size() - 1;
    } else {
        sample_index = distance(samples.begin(), sample_it);
    }

    for (uint i = 0; i != records.size(); ++i) {
        if (pos <= records[i].pos and records[i].pos + records[i].ref.length() <= pos_to and records[i].chrom == name) {
            records[i].samples[sample_index] = "0";
            //cout << "update record " << records[i] << endl;
        }
    }
    return;
}

void VCF::clear() {
    records.clear();
}

void VCF::sort_records() {
    sort(records.begin(), records.end());
    return;
}

bool VCF::pos_in_range(const uint from, const uint to) {
    for (uint i = 0; i != records.size(); ++i) {
        if (from < records[i].pos and records[i].pos + records[i].ref.length() <= to) {
            return true;
        }
    }
    return false;
}

// NB in the absence of filter flags being set to true, all results are saved. If one or more filter flags for SVTYPE are set, 
// then only those matching the filter are saved. Similarly for GRAPHTYPE.
void
VCF::save(const string &filepath, bool simple, bool complexgraph, bool toomanyalts, bool snp, bool indel, bool phsnps,
          bool complexvar) {
    /*if (samples.size() == 0)
    {
	    cout << now() << "Did not save VCF for sample" << endl;
	    return;
    }*/
    cout << now() << "Saving VCF to " << filepath << endl;

    // find date
    time_t t = time(0);
    char mbstr[10];
    strftime(mbstr, sizeof(mbstr), "%d/%m/%y", localtime(&t));

    // open and write header
    ofstream handle;
    handle.open(filepath);

    handle << "##fileformat=VCFv4.3" << endl;
    handle << "##fileDate==" << mbstr << endl;
    handle << "##ALT=<ID=SNP,Description=\"SNP\">" << endl;
    handle << "##ALT=<ID=PH_SNPs,Description=\"Phased SNPs\">" << endl;
    handle << "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">" << endl;
    handle << "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and indels\">" << endl;
    handle << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">" << endl;
    handle << "##ALT=<ID=SIMPLE,Description=\"Graph bubble is simple\">" << endl;
    handle << "##ALT=<ID=NESTED,Description=\"Variation site was a nested feature in the graph\">" << endl;
    handle
            << "##ALT=<ID=TOO_MANY_ALTS,Description=\"Variation site was a multinested feature with too many alts to include all in the VCF\">"
            << endl;
    handle << "##INFO=<ID=GRAPHTYPE,Number=1,Type=String,Description=\"Type of graph feature\">" << endl;
    handle << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (uint i = 0; i != samples.size(); ++i) {
        handle << "\t" << samples[i];
    }
    handle << endl;

    sort_records();

    for (uint i = 0; i != records.size(); ++i) {
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
    uint added = 0;
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

void VCF::write_aligned_fasta(const string &filepath, const vector<LocalNodePtr> &lmp) {
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
    uint max_len = 0;
    uint n = 0; // position in lmp
    uint ref_len = 0;
    vector<string> seqs(samples.size(), "");
    vector<uint> alt_until(samples.size(), 0);

    for (uint i = 0; i != records.size(); ++i) {
        if (records[i].pos != (uint) prev_pos) {
            // equalise lengths of sequences
            for (uint j = 0; j != samples.size(); ++j) {
                if (seqs[j].length() < max_len) {
                    //cout << "equalise length of seq " << j << endl;
                    string s(max_len - seqs[j].length(), '-'); // s == "------"
                    seqs[j] += s;
                }
            }

            // add ref sequence for gaps
            while (ref_len < records[i].pos and n < lmp.size()) {
                for (uint j = 0; j != samples.size(); ++j) {
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
        for (uint j = 0; j != samples.size(); ++j) {
            if (records[i].samples[j] == "0" and records[i].pos != (uint) prev_pos and
                pos_in_range(records[i].pos, records[i].pos + records[i].ref.length()) == false) {
                //cout << j << " add ref allele at site" << endl;
                seqs[j] += records[i].ref;
                max_len = max(max_len, (uint) seqs[j].length());
                ref_len += records[i].ref.length();
                n++;
            } else if (records[i].samples[j] == "1") {
                //cout << j << " add alt allele at site" << endl;
                assert(records[i].pos != (uint) prev_pos);
                seqs[j] += records[i].alt;
                max_len = max(max_len, (uint) seqs[j].length());
                alt_until[j] = records[i].pos + records[i].ref.length();
            }
        }
    }

    // equalise lengths of sequences
    for (uint j = 0; j != samples.size(); ++j) {
        if (seqs[j].length() < max_len) {
            //cout << "equalise length of seq " << j << endl;
            string s(max_len - seqs[j].length(), '-'); // s == "------"
            seqs[j] += s;
        }
    }

    // add ref sequence for end gaps
    while (n < lmp.size()) {
        for (uint j = 0; j != samples.size(); ++j) {
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

    for (uint j = 0; j != samples.size(); ++j) {
        handle << ">" << samples[j] << endl;
        handle << seqs[j] << endl;
    }
    handle.close();
    return;

}

bool VCF::operator==(const VCF &y) const {
    if (records.size() != y.records.size()) { return false; }
    for (uint i = 0; i != y.records.size(); ++i) {
        if (find(records.begin(), records.end(), y.records[i]) == records.end()) {
            return false;
        }
    }
    return true;
}
