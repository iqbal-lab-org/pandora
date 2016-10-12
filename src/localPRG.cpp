#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "minimizer.h"
#include "localPRG.h"
#include "interval.h"
#include "localgraph.h"
#include "errormessages.h"
//#include "path.h" // for sketching

using std::vector;
using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p, uint32_t w, uint32_t k): id(i), name(n), seq(p) 
{
    //cout << seq << endl;
    vector<uint32_t> v;
    vector<uint32_t> b = build_graph(Interval(0,seq.size()), v);
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

string LocalPRG::string_along_path(Path p)
{
    string s;
    for (deque<Interval>::iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
	s += seq.substr(it->start, it->length);
    }
    return s;
}

vector<Interval> LocalPRG::splitBySite(Interval i)
{
    // Splits interval by next_site based on substring of seq in the interval

    // Split first by var site
    vector<Interval> v;
    string::size_type k = i.start;
    string d = buff + to_string(next_site) + buff;
    string::size_type j = seq.find(d, k);
    while (j!=string::npos and j+d.size()<i.end) {
        v.push_back(Interval(k, j));
        k = j + d.size();
        j = seq.find(d, k);
    }
    //cout << "j: " << j << endl;
    if (j!=string::npos and j<i.end and j+d.size()>i.end) {
	//cout << "a1" << endl;
        v.push_back(Interval(k, j));
    } else if (j!=string::npos and j+d.size()==i.end) {
        v.push_back(Interval(k, j));
	//cout << "a2" << endl;
	if (seq.find(buff,j+d.size())==j+d.size()){
	    v.push_back(Interval(j+d.size(), j+d.size()));
	    //cout << "a2a" << endl;
	}
    } else {
	v.push_back(Interval(k, i.end));
	//cout << "a3" << endl;
    }

    // then split by var site + 1
    vector<Interval> w;
    d = buff + to_string(next_site+1) + buff;
    for(uint32_t l=0; l!=v.size(); ++l)
    {
	k = v[l].start;
	j = seq.find(d, k);
        while (j!=string::npos and j+d.size()<v[l].end) {
            w.push_back(Interval(k, j));
            k = j + d.size();
            j = seq.find(d, k);
    	}
	//cout << "j: " << j << endl;
        if (j!=string::npos and j<v[l].end and j+d.size()>v[l].end) {
	    //cout << "b1" << endl;
            w.push_back(Interval(k, j));
	} else if (j!=string::npos and j+d.size()==v[l].end) {
	    //cout << "b2" << endl;
	    w.push_back(Interval(k, j));
	    if (seq.find(buff,j+d.size())==j+d.size()){
		//cout << "b2a" << endl;
                v.push_back(Interval(j+d.size(), j+d.size()));
            }
        } else {
	    //cout << "b3" << endl;
	    w.push_back(Interval(k, v[l].end));
	}
    }
    //cout << "This was performed on sequence: " << seq << endl;
    return w;
}

vector<uint32_t> LocalPRG::build_graph(Interval i, vector<uint32_t> from_ids)
{
    // we will return the ids on the ends of any stretches of graph added
    vector<uint32_t> end_ids;

    // add nodes
    string s = seq.substr(i.start, i.length); //check length correct with this end...
    if (isalpha_string(s)) // should return true for empty string too
    {
	prg.add_node(next_id, s, i);
	// add edges from previous part of graph to start of this interval
        //cout << "from_ids: ";
        for (uint32_t j=0; j!=from_ids.size(); j++)
        {
            //cout << from_ids[j];
            prg.add_edge(from_ids[j], next_id);
        }
        //cout << endl;
	end_ids.push_back(next_id);
	next_id++;
    } else {
	// split by next var site
	vector<Interval> v = splitBySite(i); // should have length at least 4
	if(v.size()<(uint)4)
        {
	    cerr << SPLIT_SITE_WRONG_SIZE_ERROR << endl;
            cerr << "Size of partition based on site " << next_site << " is " << v.size() << endl;
            exit(-1);
        }
	next_site += 2;
	// add first interval (should be alpha)
	s = seq.substr(v[0].start, v[0].length);
	if (! (isalpha_string(s)))
        {
            cerr << SPLIT_SITE_ERROR << endl;
            cerr << "After splitting by site " << next_site << " do not have alphabetic sequence before var site: " << v[0] << endl;
            exit(-1);
        }
	prg.add_node(next_id, s, v[0]);
	// add edges from previous part of graph to start of this interval
        //cout << "from_ids: ";
        for (uint32_t j=0; j!=from_ids.size(); j++)
        {
            //cout << from_ids[j];
            prg.add_edge(from_ids[j], next_id);
        }
        //cout << endl;

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
    set<Path> walk_paths;
    Path current_path;
    for (map<uint32_t,*LocalNode>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	for (uint32_t i=it->start; i!=it->end; ++i)
	{
	    walk_paths = prg.walk(it->id, i, w+k-1);
	    for (set<Path>::iterator it2=walk_paths.begin(); it2!=walk_paths.end(); ++it)
	    {
		// find minimizer for this path	
	    }
	}
    }
}*/

/*void LocalPRG::minimizer_sketch (uint32_t w, uint32_t k)
{
    set<Path> current_paths, new_paths, paths;

    // if have no nodes, just return, can't get minimizers
    if (prg.nodes.size()<=0)
	return;

    current_paths.insert(Path({prg.nodes.at(0)->pos});
    while (current_paths.size()>0)
    {
        for (std::set<Path>::iterator it=current_paths.begin(); it!=current_paths.end(); ++it)
	{
	    if (it->length < w+k-1 + it->path.begin()->length - 1) // i.e. if any of the paths 
	    {
		// if too short increase path size
		paths = prg.extend_path(*it);
		new_paths.insert(paths.begin(), paths.end());
	    } else {
	        //find minimizers for all paths starting in first interval of the path
	        for (uint32_t j=it->path.begin()->start; j!=it->path.begin()->end; ++j)
		{
		}
	    }

	}
    }
}*/

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

