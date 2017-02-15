#include <iostream>
#include <cassert>
#include <cmath>
#include "maxpath.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

MaxPath::MaxPath(){};

MaxPath::MaxPath(vector<LocalNode*> x, vector<int> y, uint32_t z): kmers_on_path(y), num_equivalent_paths(z)
{
    npath.reserve(100);
    npath = x;
}

void MaxPath::extend(const MaxPath new_mp)
{
    uint old_size = npath.size();
    vector<LocalNode*>::const_iterator it = new_mp.npath.begin();
    while((*it)->id <= npath.back()->id and it!=new_mp.npath.end())
    {
	it++;
    }
    npath.insert(npath.end(),it,new_mp.npath.end());
    assert(npath.size()>old_size);

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
    // instead return mean log prob (not log mean prob)
    float p = 0;
    uint t = 0;
    assert(kmers_on_path.size() == kmer_path_probs.size());

    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == 1)
        {
            p += kmer_path_probs[i];
            t += 1;
        }
    }
    if (t == 0)
    {
	mean_prob = 0;
	return 1;
    }
    mean_prob = p/t;
    return p/t;
}

bool MaxPath::has_at_least_n_hit_minis_on_path(const vector<uint32_t>& counts, uint32_t n)
{
    assert(kmers_on_path.size() == counts.size() || assert_msg("kmers_on_path.size(): " << kmers_on_path.size() << ", counts.size(): " << counts.size()));

    uint tally = 0;
    for (uint i = 0; i!=kmers_on_path.size(); ++i)
    {
        if (kmers_on_path[i] == 1 and counts[i] >= 1)
        {
	    tally++;
	    if (tally >= n)
            {
		return true;
	    }
        }
    }
    return false;
}

bool VMPgreater::operator()( const vector<MaxPath>& lx, const vector<MaxPath>& rx )
{
    return lx[0].prob > rx[0].prob;
};
