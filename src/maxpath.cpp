#include <iostream>
#include <cassert>
#include "maxpath.h"

using namespace std;

MaxPath::MaxPath(){};

//MaxPath::MaxPath(vector<LocalNode*> x, vector<bool> y, uint32_t z): npath(x), kmers_on_path(y), num_equivalent_paths(z){};
MaxPath::MaxPath(vector<LocalNode*> x, vector<bool> y, uint32_t z)
{
    cout << "make maxpath" << endl;
    npath = x;
    cout << "added npath" << endl;
    kmers_on_path = y;
    cout << "added kmers_on_path" << endl;
    num_equivalent_paths = z;
    cout << "done creating maxpath" << endl;
}
    

void MaxPath::extend(const MaxPath new_mp)
{
    uint old_size = npath.size();
    npath.insert(npath.end(),new_mp.npath.begin(),new_mp.npath.end());
    assert(npath.size() == old_size + new_mp.npath.size());
    assert(new_mp.kmers_on_path.size()==kmers_on_path.size());
    for (uint i=0; i!=new_mp.kmers_on_path.size(); ++i)
    {
	if (new_mp.kmers_on_path[i]==true)
	{
	    kmers_on_path[i] = true;
	}
    }
    return;
}

float MaxPath::get_prob(const vector<float>& kmer_path_probs)
{
    float p = 0;
    assert(kmers_on_path.size() == kmer_path_probs.size());

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == true)
        {
            p += kmer_path_probs[i];
        }
    }
    prob = p;
    return p;
}

float MaxPath::get_mean_prob(const vector<float>& kmer_path_probs)
{
    float p = 0;
    uint t = 0;
    assert(kmers_on_path.size() == kmer_path_probs.size());

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == true)
        {
            p += kmer_path_probs[i];
            t += 1;
        }
    }
    if (t == 0)
    {
	mean_prob = 0;
	return 0;
    }
    mean_prob = p/t;
    return p/t;
}

