#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include "inthash.h"
#include "minimizer.h"
#include "seq.h"
#include "utils.h"

using std::vector;
using namespace std;


Seq::Seq (uint32_t i, string n, string p, uint32_t w, uint32_t k): id(i), name(n), seq(p) {
    minimizer_sketch (w, k);
}

Seq::~Seq()
{
    for (auto c : sketch)
    {
        delete c;
    }
}

void Seq::minimizer_sketch (uint32_t w, uint32_t k)
{
    // If sequence too short, just return
    if (seq.length()+1 < w+k) {//cout << "Sequence too short to sketch" << endl; 
	return;}

    // Store minimizers in vector until the end, so don't have to sort every time add a new element
    deque<Minimizer*> v;

    // for each window position
    for(uint32_t wpos=0; wpos <= seq.length()-w-k+1 ; ++wpos)
    {
	// find the lex smallest kmer in the window
	//set<string> kmers;
	set<uint64_t> kmers;
	string kmer;
	uint64_t kh;

	for (uint32_t i = 0; i < w; i++)
	{
	    kmer = seq.substr(wpos+i, k);
	    kh = kmerhash(kmer, k);
	    kmers.insert(kh);
	}
	uint64_t smallest_word = *kmers.begin();
	for (uint32_t i = 0; i < w; i++)
	{
	    kmer = seq.substr(wpos+i, k);
	    kh = kmerhash(kmer, k);
	    if (kh == smallest_word)
            {
            	Minimizer *m;
		m = new Minimizer(kh, wpos+i, wpos+i+k);
 		pointer_values_equal<Minimizer> eq = { m };
                if ( find_if(v.begin(), v.end(), eq) == v.end() )
		{
		    v.push_front(m);
		} else {
		    delete m;
		}
	    }
	}
    }
    
    // Now add the vector of minimizers to the sketch set
    sketch.insert(v.begin(), v.end());
    return;
}

std::ostream& operator<< (std::ostream & out, Seq const& data) {
    out << data.name;
    return out ;
}

