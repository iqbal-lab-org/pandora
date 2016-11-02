#include <iostream>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include "minimizer.h"
#include "localPRG.h"
#include "interval.h"
#include "localgraph.h"
#include "index.h"
#include "errormessages.h"
#include "path.h"
#include "minihits.h"
#include "minihit.h"
#include "inthash.h"

using std::vector;
using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p): next_id(0),buff(" "), next_site(5), id(i), name(n), seq(p)
{
    //cout << seq << endl;
    vector<uint32_t> v;
    // avoid error if a prg contains only empty space as it's sequence
    if(seq.find_first_not_of("\t\n\v\f\r") != std::string::npos)
    {
        vector<uint32_t> b = build_graph(Interval(0,seq.size()), v);
        //minimizer_sketch (idx, w, k);
    } else {
	prg.add_node(0, "", Interval(0,0));
    }
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
	    //cout << "Found non-alpha char: " << s[j] << endl;
            return 0; //False
        }
    return 1; //True
}

string LocalPRG::string_along_path(Path p)
{
    //cout << p << endl;
    assert(p.start<=seq.length());
    assert(p.end<=seq.length());
    string s;
    for (deque<Interval>::iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
	s += seq.substr(it->start, it->length);
        //cout << s << endl;
    }
    //cout << "lengths: " << s.length() << ", " << p.length << endl;
    assert(s.length()==p.length);
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

void LocalPRG::minimizer_sketch (Index* idx, uint32_t w, uint32_t k)
{
    //cout << "START SKETCH FUNCTION" << endl;
    vector<Path> walk_paths;
    Path current_path, kmer_path;
    string kmer;
    uint64_t kh, smallest = std::numeric_limits<uint64_t>::max();;
    for (map<uint32_t,LocalNode*>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	for (uint32_t i=it->second->pos.start; i!=it->second->pos.end; ++i)
	{
	    walk_paths = prg.walk(it->second->id, i, w+k-1);
            //cout << "for id, i: " << it->second->id << ", " << i << " found " << walk_paths.size() << " paths" << endl;
	    for (vector<Path>::iterator it2=walk_paths.begin(); it2!=walk_paths.end(); ++it2)
	    {
		//cout << "Minimize path: " << *it2 << endl;
		// find minimizer for this path	
		smallest = std::numeric_limits<uint64_t>::max();
		for (uint32_t j = 0; j != w; j++)
		{
		    //cout << "i+j" << i+j << endl;
		    kmer_path = it2->subpath(i+j,k);
		    if (kmer_path.path.size() > 0)
		    {
                        //cout << "found path" << endl;
		        kmer = string_along_path(kmer_path);
                        kh = kmerhash(kmer, k);
		        smallest = min(smallest, kh);
		    }
		}
                //cout << "smallest word: " << smallest_word << endl;
		for (uint32_t j = 0; j != w; j++)
                {
                    //cout << "i+j" << i+j << endl;
                    kmer_path = it2->subpath(i+j,k);
		    //cout << "kmer path" << kmer_path << endl;
		    if (kmer_path.path.size() > 0)
		    {
                    	//cout << "found path" << endl;
                    	kmer = string_along_path(kmer_path);
			kh = kmerhash(kmer, k);
		    	if (kh == smallest)
		    	{
                            //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
			    idx->add_record(kh, id, kmer_path);
		    	}
		    }
                }
	    }
	}
    }
}

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

void LocalPRG::get_covgs(MinimizerHits* minimizer_hits)
{
    //first, for each localnode of the localgraph, set covg to 0
    for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
    {
	n->second->covg = 0;
    }
    // for each hit
    for (set<MinimizerHit*, pComp>::iterator mh = minimizer_hits->hits.begin(); mh != minimizer_hits->hits.end(); ++mh)
    {
        // if it is against this prg id
        if ((*mh)->prg_id == id)
        {
            // then for each interval of the hit
            for (deque<Interval>::const_iterator it=(*mh)->prg_path.path.begin(); it!=(*mh)->prg_path.path.end(); ++it)
            {
                // update the covg on the appropriate node of the prg
                for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
		{
                    if (it->end > n->second->pos.start and it->start < n->second->pos.end)
            	    {
                        n->second->covg += min(it->end, n->second->pos.end) - max(it->start, n->second->pos.start);
                    } else if (it->end < n->second->pos.start) 
		    { break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
		}
            }
        }
    }
    return;
}

void LocalPRG::update_covg_with_hit(MinimizerHit* mh)
{
    // then for each interval of the hit
    for (deque<Interval>::const_iterator it=mh->prg_path.path.begin(); it!=mh->prg_path.path.end(); ++it)
    {
        // update the covg on the appropriate node of the prg
        for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
        {
            if (it->end > n->second->pos.start and it->start < n->second->pos.end)
            {
                n->second->covg += min(it->end, n->second->pos.end) - max(it->start, n->second->pos.start);
            } else if (it->end < n->second->pos.start)
            { break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
        }
    }
}
