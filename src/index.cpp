#include <iostream>
#include <fstream>
#include <cassert>
#include <map>
#include <vector>
#include "minirecord.h"
#include "path.h"
#include "index.h"
#include "utils.h"

using namespace std;

Index::Index() {};

Index::~Index() {
    clear();
};

void Index::add_record(uint64_t kmer, uint32_t prg_id, Path path, bool strand)
{
    //cout << "Add kmer " << kmer << " id, path " << prg_id << ", " << path << endl;
    map<uint64_t, vector<MiniRecord>>::iterator it=minhash.find(kmer);
    if(it==minhash.end())
    {
        vector<MiniRecord> newv;
        newv.reserve(20);
	newv.push_back(MiniRecord(prg_id, path, strand));
        minhash.insert(pair<uint64_t, vector<MiniRecord>>(kmer,newv));
        //cout << "New minhash size: " << minhash.size() << endl; 
    } else {
	MiniRecord mr(prg_id, path, strand);	
        if (find(it->second.begin(), it->second.end(), mr)==it->second.end())
	{
            it->second.push_back(mr);
	}
        //cout << "New minhash entry for  kmer " << kmer << endl;
    }
}

void Index::clear()
{
    for(auto it = minhash.begin(); it != minhash.end();)
        it = minhash.erase(it);
}

void Index::save(const string& prgfile)
{
    cout << now() << "Saving index" << endl;
    ofstream handle;
    handle.open (prgfile + ".idx");

    for (auto it = minhash.begin(); it != minhash.end(); ++it)
    {
        handle << it->first;
        for (uint j = 0; j!=it->second.size(); ++j)
        {
            handle << "\t" << it->second[j];
        }
        handle << endl;

    }
    handle.close();
    cout << now() << "Finished saving " << minhash.size() << " entries to file" << endl;
    return;
}

void Index::load(const string& prgfile)
{
    cout << now() << "Loading index" << endl;
    //string line;
    //vector<string> vstring;
    uint32_t key;
    MiniRecord mr;
    vector<MiniRecord> vmr;

    ifstream myfile (prgfile + ".idx");
    if (myfile.is_open())
    {
	while (myfile.good())
	{
	    if (myfile >> key)
	    {
		minhash[key] = vmr;
	    } else if (myfile >> mr) {
	        minhash[key].push_back(mr);
	    }
	}
        /*while ( getline (myfile,line).good() )
	{
	    vstring = split(line, "\t");
	    assert(vstring.size() >= 2);
	    vstring[1] << >> mr;
	    minhash[(uint)vstring[0]] = {mr};
	    for (uint i=2; i!=vstring.size(); ++i)
	    {
		vstring[i] << >> mr;
		minhash[(uint)vstring[0]].push_back(mr);
	    }
	}
	//last line
	vstring = split(line, "\t");
        assert(vstring.size() >= 2);
        vstring[1] << >> mr;
        minhash[(uint)vstring[0]] = {mr};
        for (uint i=2; i!=vstring.size(); ++i)
        {
            vstring[i] << >> mr;
            minhash[(uint)vstring[0]].push_back(mr);
        }*/
    } else {
        cerr << "Unable to open index file " << prgfile << ".idx" << endl;
	exit(1);
    }
    cout << now() << "Finished loading " << minhash.size() << " entries to index" << endl;
    return;
}
	    

