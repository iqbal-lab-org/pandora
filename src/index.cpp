#include <iostream>
#include <map>
#include <vector>
#include "minirecord.h"
#include "path.h"
#include "index.h"

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
