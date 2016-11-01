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
    //cout << "Start sketching" << endl;
    // If sequence too short, just return
    if (seq.length()+1 < w+k) {//cout << "Sequence too short to sketch" << endl; 
	return;}

    // initializations
    string kmer;
    uint64_t kh, smallest;
    Minimizer *m;
    Minimizer *m_previous;

    // for each window position
    for(uint32_t wpos=0; wpos <= seq.length()-w-k+1 ; ++wpos)
    {
	//cout << "wpos: " << wpos << endl;
    	// if wpos==0 or the previous minimizer starts outside current window, calculate new mini from scratch
        if((wpos == 0) or (m_previous->pos.start < wpos))
	{
            smallest = std::numeric_limits<uint64_t>::max();
	    // find the lex smallest kmer in the window
	    for (uint32_t i = 0; i < w; i++)
	    {
	    	kmer = seq.substr(wpos+i, k);
	        kh = kmerhash(kmer, k);
	        smallest = min(smallest, kh);
	    }
	    for (uint32_t i = 0; i < w; i++)
	    {
	        kmer = seq.substr(wpos+i, k);
	        kh = kmerhash(kmer, k);
	        if (kh == smallest)
                {
		    m = new Minimizer(kh, wpos+i, wpos+i+k);
		    sketch.insert(m);
		    m_previous = m;
	        }
	    }
        } else {
        // otherwise only need to do something if the kmer from the newest position at end of new window is smaller or equal to the previous smallest
	    kmer = seq.substr(wpos+w-1, k);
            kh = kmerhash(kmer, k);
	    //cout << "Last kh for wpos: " << kh << " compared to previous smallest: " << smallest << endl;
	    if(kh <= smallest)
	    {
	        m = new Minimizer(kh, wpos+w-1, wpos+w-1+k);
		sketch.insert(m);
		m_previous = m;
       	    }
            smallest = min(smallest, kh);
	}
    }
    
    return;
}

std::ostream& operator<< (std::ostream & out, Seq const& data) {
    out << data.name;
    return out ;
}

