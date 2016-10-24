#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
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


    // for each window position
    for(uint32_t wpos=0; wpos <= seq.length()-w-k+1 ; ++wpos)
    {
	// find the lex smallest kmer in the window
	set<string> kmers;
	string kmer;

	for (uint32_t i = 0; i < w; i++)
	{
	    kmer = seq.substr(wpos+i, k);
	    kmers.insert(kmer);
	}
	string smallest_word = *kmers.begin();
	for (uint32_t i = 0; i < w; i++)
	{
	    kmer = seq.substr(wpos+i, k);
	    if (kmer == smallest_word)
            {
            	Minimizer *m;
		m = new Minimizer(kmer, wpos+i, wpos+i+k);
 		pointer_values_equal<Minimizer> eq = { m };
                if ( find_if(sketch.begin(), sketch.end(), eq) == sketch.end() )
		{
		    sketch.insert(m);
		} else {
		    delete m;
		}
	    }
	}
    }
    return;
}

std::ostream& operator<< (std::ostream & out, Seq const& data) {
    out << data.name;
    return out ;
}

