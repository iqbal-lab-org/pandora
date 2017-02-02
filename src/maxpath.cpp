#include <iostream>
#include <cassert>
#include <cmath>
#include "maxpath.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

MaxPath::MaxPath(){};

MaxPath::MaxPath(vector<LocalNode*> x, vector<int> y, uint32_t z): npath(x), kmers_on_path(y), num_equivalent_paths(z)
{
    npath.reserve(100);
}
/*MaxPath::MaxPath(vector<LocalNode*> x, vector<bool> y, uint32_t z)
{
    cout << "make maxpath" << endl;
    npath = x;
    cout << "added npath" << endl;
    kmers_on_path = y;
    cout << "added kmers_on_path" << endl;
    num_equivalent_paths = z;
    cout << "done creating maxpath" << endl;
}*/
    

void MaxPath::extend(const MaxPath new_mp)
{
    cout << "extend npaths" << endl;
    uint old_size = npath.size();
    cout << "old size: " << old_size << endl;
    vector<LocalNode*>::const_iterator it = new_mp.npath.begin();
    while((*it)->id <= npath.back()->id and it!=new_mp.npath.end())
    {
	it++;
    }
    /*for (uint n=0; n!=npath.size(); ++n)
    {
	if (*npath[n] == **it)
	{
	    it++;
	}
    }*/
    npath.insert(npath.end(),it,new_mp.npath.end());
    assert(npath.size()>old_size);
    cout << "new size: " << npath.size() << endl;

    cout << "take intersection of kmers" << endl;
    assert(new_mp.kmers_on_path.size()==kmers_on_path.size());
    // keep the intersection of kmers
    for (uint n=0; n!=kmers_on_path.size(); ++n)
    {
	if (new_mp.kmers_on_path[n] != kmers_on_path[n])
	{
	    kmers_on_path[n] = 0;
	}
    }
    return;
}

float MaxPath::get_prob(const vector<float>& kmer_path_probs)
{
    float p = 0;
    assert(kmers_on_path.size() == kmer_path_probs.size() || assert_msg("kmers_on_path.size(): " << kmers_on_path.size() << ", kmer_path_probs.size(): " << kmer_path_probs.size()));

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == 1)
        {
            p += kmer_path_probs[i];
        }
    }
    prob = p;
    return p;
}

float MaxPath::get_mean_prob(const vector<float>& kmer_path_probs)
{
    // currently does not give mean, testing
    float p = 0;
    uint t = 0;
    assert(kmers_on_path.size() == kmer_path_probs.size());

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == 1)
        {
            p += exp(kmer_path_probs[i]);
            t += 1;
        }
    }
    if (t == 0)
    {
	mean_prob = 0;
	return 0;
    }
    mean_prob = log(p/t);
    return log(p/t);
}

float MaxPath::get_median_prob(const vector<float>& kmer_path_probs)
{
    vector<float> ps;
    assert(kmers_on_path.size() == kmer_path_probs.size());

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == 1)
        {
            ps.push_back(kmer_path_probs[i]);
        }
    }
    if (ps.size() == 0)
    {
	median_prob = 0;
	return 0;
    }
    uint t = (ps.size() + 1)/2;
    median_prob = ps[t];
    return ps[t];
}

bool VMPgreater::operator()( const vector<MaxPath>& lx, const vector<MaxPath>& rx )
{
    return lx[0].prob > rx[0].prob;
};
