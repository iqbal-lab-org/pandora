#include <iostream>
#include <cassert>
//#include <cstring>
#include <map>
#include <vector>
#include "minirecord.h"
#include "path.h"
#include "index.h"

using namespace std;

Index::Index() {};

Index::~Index() {
    clear();
    /*for (auto c: minhash)
    {
        delete c.second;
    }*/
};

void Index::add_record(uint64_t kmer, uint32_t prg_id, Path path)
{
    //cout << "Add kmer " << kmer << " id, path " << prg_id << ", " << path << endl;
    //assert(kmer.length()==path.length);
    vector<MiniRecord> newv;
    map<uint64_t, vector<MiniRecord>>::iterator it=minhash.find(kmer);
    if(it==minhash.end())
    {
        newv.clear();
        MiniRecord mr(prg_id, path);
        newv.push_back(mr);
        minhash.insert(pair<uint64_t, vector<MiniRecord>>(kmer,newv));
        //cout << "New minhash size: " << minhash.size() << endl; 
    } else {
	MiniRecord mr(prg_id, path);	
        if (find(it->second.begin(), it->second.end(), mr)==it->second.end())
	{
            it->second.push_back(mr);
	}
        //cout << "New minhash[" << kmer << "] size: " << minhash[kmer].size() << endl;
    }
}

void Index::clear()
{
    for(auto it = minhash.begin(); it != minhash.end();)
        it = minhash.erase(it);
}
