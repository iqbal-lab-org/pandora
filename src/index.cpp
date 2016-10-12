#include <iostream>
#include <cassert>
#include <cstring>
#include <map>
#include <vector>
#include "minirecord.h"
#include "path.h"
#include "index.h"

using namespace std;

Index::Index() {};

Index::~Index() {
    /*for (auto c: minhash)
    {
        delete c.second;
    }*/
};

void Index::add_record(string miniKmer, uint32_t prg_id, Path path)
{
    assert(miniKmer.length()==path.length);
    vector<MiniRecord> newv;
    map<string, vector<MiniRecord>>::iterator it=minhash.find(miniKmer);
    if(it==minhash.end())
    {
        newv.clear();
        newv.push_back(MiniRecord(prg_id, path));
        minhash.insert(pair<string, vector<MiniRecord>>(miniKmer,newv));
    } else {
	MiniRecord mr = MiniRecord(prg_id, path);	
        if (find(it->second.begin(), it->second.end(), mr)==it->second.end())
	{
            it->second.push_back(mr);
	}
    }
}
