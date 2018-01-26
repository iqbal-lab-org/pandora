#include <iostream>
#include <cassert>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include "inthash.h"
#include "minimizer.h"
#include "seq.h"
#include "utils.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

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

void Seq::initialize(uint32_t i, string n, string p, uint32_t w, uint32_t k)
{
    id = i;
    name = n;
    seq = p;
    for (auto c : sketch)
    {   
        delete c;
    }
    sketch.clear();
    minimizer_sketch (w, k);
}

void Seq::minimizer_sketch (const uint32_t w, const uint32_t k)
{
    if (seq.length()+1 < w+k) {//cout << "Sequence too short to sketch" << endl; 
        return;}

    // initializations
    uint c, i, i_, num_minis_found=0;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, smallest, kmer[2] = {0,0}, kh[2] = {0,0};
    vector<Minimizer*> vm;
    vm.reserve(w); //can't reserve a deque
    Minimizer* m;

    char myArray[seq.size()+1];//as 1 char space for null is also required
    strcpy(myArray, seq.c_str());

    for(uint32_t buff=0; buff < seq.length() ; ++buff)
    {
	c = nt4((uint8_t)seq[buff]);
        if (c < 4) { // not an ambiguous base
            kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
            kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer	
	    kh[0] = hash64(kmer[0], mask);
            kh[1] = hash64(kmer[1], mask);
	} else {
	    cout << now() << "bad letter - found a non AGCT base in read so skipping read " << name << endl;
	    for (auto c : sketch)
    	    {   
        	delete c;
    	    }
	    sketch.clear();
	    num_minis_found = 0;
	    break;
        }
	
	if (buff >=k-1)
	{
	    //make a mini
	    m = new Minimizer(min(kh[0], kh[1]), buff-k+1, buff+1, (kh[0]<=kh[1]));
	    vm.push_back(m);
	    //cout << "kmer pos " << buff-k+1 << " and vm.size()==" << vm.size() << endl;
	}

	if (vm.size()==w)
	{
	    smallest = std::numeric_limits<uint64_t>::max();
	    // find smallest khash value for the w kmers
	    //cout << "find smallest of kmers starting at pos ";
	    for (i=0; i!=vm.size(); ++i)
	    {	
		if (vm[i]->kmer <= smallest)
		{
		    smallest = vm[i]->kmer;
		    i_ = i;
		}
		//cout << vm[i]->pos.start << " ";
	    }
	    //cout << endl;

	    // add these minimizers to sketch, and delete others
	    for (i=0; i!=vm.size(); ++i)
            {
		if (vm[i]->kmer == smallest)
		{
		    sketch.insert(vm[i]);
		    //cout << "add minikmer for pos: " << vm[i]->pos.start << " " << *vm[i] << " as new smallest, since vm.size()==" << vm.size() << endl;
		    num_minis_found += 1;
		} else if (i < i_) 
		{
		    delete vm[i];
		}
	    }

	    vm.erase(vm.begin(), vm.begin()+i_+1);
            /*for (i=0; i!=i_+1; ++i)
            {
		vm.pop_front();
	    }*/
	    
	    assert(vm.size() < w || assert_msg("we can't have added a smallest kmer correctly as vm still has size " << vm.size()));
	} else if ( buff >= w+k-1 and vm.back()->kmer <= smallest)
	{
	    sketch.insert(vm.back());
	    //cout << "add minikmer for pos: " << vm.back()->pos.start << " " << *vm.back() << " which is smaller than " << smallest << endl;
	    num_minis_found += 1;
	    smallest = vm.back()->kmer;
	    // delete all but the last kmer
	    for (i=0; i!=vm.size()-1; ++i)
	    {
		delete vm[i];
	    }
	    vm.clear();
	}

	if (buff == seq.length()-1)
        {
	    // delete remaining elements of vm
	    //cout << "we have reached final position " << buff << " and vm.size()==" << vm.size() << endl;
	    for (i=0; i!=vm.size(); ++i)
            {
                delete vm[i];
            }
            vm.clear();
	}
    }
    //cout << "for read " << name << " num_minis_found " << num_minis_found << " and sketch size " << sketch.size() << endl;
    assert(num_minis_found == sketch.size());
    //cout << now() << "Sketch size " << sketch.size() << " for read " << name << endl;
    return;
}

std::ostream& operator<< (std::ostream & out, Seq const& data) {
    out << data.name;
    return out ;
}

