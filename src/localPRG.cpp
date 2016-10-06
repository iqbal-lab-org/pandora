#include <iostream>
#include <string>
#include <vector>
//#include <set>
#include "minimizer.h"
#include "localPRG.h"
#include "interval.h"
#include "localgraph.h"
//#include "path.h" // for sketching

using std::vector;
using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p, uint32_t w, uint32_t k): id(i), name(n), seq(p) 
{
    vector<uint32_t> v;
    build_graph(Interval(0,seq.size()), v);
    //minimizer_sketch (w, k);
}

LocalPRG::~LocalPRG()
{
    for (auto c : sketch)
    {
        delete c;
    }
}

bool LocalPRG::isalpha_string ( string s )
{
    // Returns if a string s is entirely alphabetic
    for(uint32_t j=0; j<s.length(); ++j)
        if(isalpha (s[j]) == 0)
        {
            return 0; //False
        }
    return 1; //True
}

vector<Interval> LocalPRG::splitBySite(Interval i)
{
    // Splits interval by next_site based on substring of seq in the interval

    // Split first by var site
    vector<Interval> v;
    string::size_type k = i.start;
    string d = buff + to_string(next_site) + buff;
    string::size_type j = seq.find(d, k);

    if (j == string::npos or j>=i.end) {
        v.push_back(i);
    }
    while (j != string::npos and j<i.end) {
        v.push_back(Interval(k, j-k+1));
        k = j + d.size();
        j = seq.find(d, k);

        if (j == string::npos or j>=i.end)
            v.push_back(Interval(k, i.end));
    }

    // then split by var site + 1
    vector<Interval> w;
    d = buff + to_string(next_site+1) + buff;
    for(uint32_t l=0; l<v.size(); ++l)
    {
	k = v[l].start;
	j = seq.find(d, k);
	
	if (j == string::npos or j>=v[l].end) {
            w.push_back(v[l]);
    	}
        while (j != string::npos and j<v[l].end) {
            w.push_back(Interval(k, j-k+1));
            k = j + d.size();
            j = seq.find(d, k);

            if (j == string::npos or j>=v[l].end)
            	w.push_back(Interval(k, v[l].end));
    	}
    }
    return w;
}

vector<uint32_t> LocalPRG::build_graph(Interval i, vector<uint32_t> from_ids)
{
    // we will return the ids on the ends of any stretches of graph added
    vector<uint32_t> end_ids;

    // add edges from previous part of graph to start of this interval
    for (uint32_t j=0; j!=from_ids.size(); j++)
    {
        prg.add_edge(from_ids[j], next_id);
    }

    // add nodes
    string s = seq.substr(i.start, i.length); //check length correct with this end...
    if (isalpha_string(s)) // should return true for empty string too
    {
	prg.add_node(next_id, s, i);
	end_ids.push_back(next_id);
	next_id++;
    } else {
	// split by next var site
	vector<Interval> v = splitBySite(i); // should have length at least 4
	next_site += 2;
	// add first interval (should be alpha)
	s = seq.substr(v[0].start, v[0].length);
	prg.add_node(next_id, s, v[0]);
	vector<uint32_t> mid_ids;
	mid_ids.push_back(next_id);
	next_id++;
	// add (recurring as necessary) middle intervals
	for (uint32_t j=1; j!=v.size() - 1; j++)
	{
	    vector<uint32_t> w = build_graph(v[j], mid_ids);
	    end_ids.insert(end_ids.end(), w.begin(), w.end());
	}
	// add end interval
	end_ids = build_graph(v.back(), end_ids);
    }
    return end_ids;
}

/*void LocalPRG::minimizer_sketch (uint32_t w, uint32_t k)
{
    // If sequence too short, just return
    if (seq.length()+1 < w+k) {cout << "Sequence too short to sketch" << endl; return;}

    // for each window position
    for(uint32_t wpos=0; wpos <= seq.length()-w-k+1 ; ++wpos)
    {
	// find the lex smallest kmer in the window
	set<string> kmers;

	for (uint32_t i = 0; i < w; i++)
	{
	    string kmer = seq.substr(wpos+i, k);
	    kmers.insert(kmer);
	}
	string smallest_word = *kmers.begin();
	for (uint32_t i = 0; i < w; i++)
	{
	    string kmer = seq.substr(wpos+i, k);
	    if (kmer == smallest_word)
            {
		deque<interval> d = {interval(wpos+i, wpos+i+k)};
            	Minimizer *m;
		m = new Minimizer(kmer, d);
		sketch.push_back(m);
	    }
	}
    }

    //cout << "Found " << sketch.size() << " minimizers." << endl;
    sort(sketch.begin(), sketch.end(), pMiniComp);
    return;
}*/

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

