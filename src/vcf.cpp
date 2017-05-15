#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include "vcfrecord.h"
#include "vcf.h"
#include "utils.h"

using namespace std;

VCF::VCF() {};

VCF::~VCF() {
    clear();
};

void VCF::add_record(std::string c, uint32_t p, std::string r, std::string a)
{
    VCFRecord vr(c, p, r, a);
    if ( find(records.begin(), records.end(), vr) == records.end())
    {
	records.push_back(vr);
    }
}

void VCF::add_record(VCFRecord& vr)
{
    if ( find(records.begin(), records.end(), vr) == records.end())
    {
        records.push_back(vr);
    }
}

void VCF::clear()
{
    records.clear();
}

void VCF::save(const string& filepath)
{
    cout << now() << "Saving VCF to " << filepath << endl;
    ofstream handle;
    handle.open (filepath);

    handle << "##fileformat=VCFv4.3" << endl;
    handle << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" << endl;

    sort(records.begin(), records.end()); // we need the records in order for it to be a valid vcf

    for (uint i=0; i!=records.size(); ++i)
    {
        handle << records[i] << endl;
    }
    handle.close();
    cout << now() << "Finished saving " << records.size() << " entries to file" << endl;
    return;
}

void VCF::load(const string& filepath)
{
    cout << now() << "Loading VCF from " << filepath << endl;
    VCFRecord vr;
    string line;
    stringstream ss;
    uint added = 0;
    // NB this doesn't currently clear records first. Do we want to?

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            if (line[0] != '#')
            {
		ss << line;
                ss >> vr;
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

