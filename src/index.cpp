#include <iostream>
#include <fstream>
#include <cassert>
#include <unordered_map>
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
    //cout << "Add kmer " << kmer << " id, path, strand " << prg_id << ", " << path << ", " << strand << endl;
    unordered_map<uint64_t, vector<MiniRecord>*>::iterator it=minhash.find(kmer);
    if(it==minhash.end())
    {
        vector<MiniRecord>* newv = new vector<MiniRecord>;
        newv->reserve(20);
	newv->push_back(MiniRecord(prg_id, path, strand));
        minhash.insert(pair<uint64_t, vector<MiniRecord>*>(kmer,newv));
        //cout << "New minhash size: " << minhash.size() << endl; 
    } else {
	MiniRecord mr(prg_id, path, strand);	
        if (find(it->second->begin(), it->second->end(), mr)==it->second->end())
	{
            it->second->push_back(mr);
	}
        //cout << "New minhash entry for  kmer " << kmer << endl;
    }
}

void Index::clear()
{
    for(auto it = minhash.begin(); it != minhash.end();)
    {
	delete it->second;
        it = minhash.erase(it);
    }
}

void Index::save(const string& prgfile, uint32_t w, uint32_t k)
{
    cout << now() << "Saving index" << endl;
    ofstream handle;
    handle.open (prgfile + ".k" + to_string(k) + ".w" + to_string(w) + ".idx");

    handle << minhash.size() << endl;

    for (auto it = minhash.begin(); it != minhash.end(); ++it)
    {
        handle << it->first << "\t" << it->second->size();
        for (uint j = 0; j!=it->second->size(); ++j)
        {
            handle << "\t" << (*(it->second))[j];
        }
        handle << endl;

    }
    handle.close();
    cout << now() << "Finished saving " << minhash.size() << " entries to file" << endl;
    return;
}

void Index::load(const string& prgfile, uint32_t w, uint32_t k)
{
    cout << now() << "Loading index" << endl;
    //string line;
    //vector<string> vstring;
    uint32_t key;
    size_t size;
    int c;
    MiniRecord mr;
    bool first = true;
    //vector<MiniRecord> vmr;

    ifstream myfile (prgfile + ".k" + to_string(k) + ".w" + to_string(w) + ".idx");
    if (myfile.is_open())
    {
	while (myfile.good())
	{
	    c = myfile.peek();
	    if (isdigit(c) and first == true)
	    {
		myfile >> size;
		minhash.reserve(size);
		first = false;
		myfile.ignore(1,'\t');
	    } else if (isdigit(c) and first == false) {
		myfile >> key;
		myfile.ignore(1,'\t');
		myfile >> size;
		vector<MiniRecord>* vmr = new vector<MiniRecord>;
		vmr->reserve(size);
		minhash[key] = vmr;
		myfile.ignore(1,'\t');
	    } else {
		myfile >> mr;
	        minhash[key]->push_back(mr);
		myfile.ignore(1,'\t');
	    }
	}
    } else {
        cerr << "Unable to open index file " << prgfile << ".idx" << endl;
	exit(1);
    }
    cout << now() << "Finished loading " << minhash.size() << " entries to index" << endl;
    return;
}
	    

