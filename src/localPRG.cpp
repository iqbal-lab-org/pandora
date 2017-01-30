#include <iostream>
#include <fstream>
#include <algorithm>
#include <stdint.h>
#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <cassert>
#include <limits>
#include <cmath>
#include "minimizer.h"
#include "localPRG.h"
#include "interval.h"
#include "localgraph.h"
#include "index.h"
#include "path.h"
#include "minihits.h"
#include "minihit.h"
#include "inthash.h"
#include "pannode.h"
#include "utils.h"
#include "maxpath.h"
#include "errormessages.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using std::vector;
using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p): next_id(0),buff(" "), next_site(5), id(i), name(n), seq(p)//, max_level(0), num_minis(0)
{
    //cout << "Making new LocalPRG instance " << name << endl;
    vector<uint32_t> v;
    // avoid error if a prg contains only empty space as it's sequence
    if(seq.find_first_not_of("\t\n\v\f\r") != std::string::npos)
    {
	//cout << "build PRG graph ...";
        vector<uint32_t> b = build_graph(Interval(0,seq.size()), v);
	//cout << " done! Found " << prg.nodes.size() << " nodes" << endl;
        //minimizer_sketch (idx, w, k);
    } else {
	//cout << "seq empty" << endl;
	prg.add_node(0, "", Interval(0,0));
    }
}

/*LocalPRG::~LocalPRG()
{
    for (auto c : sketch)
    {
        delete c;
    }
}*/

bool LocalPRG::isalpha_string (const string& s )
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

string LocalPRG::string_along_path(const Path& p)
{
    //cout << p << endl;
    assert(p.start<=seq.length());
    assert(p.end<=seq.length());
    string s;
    for (deque<Interval>::const_iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
	s += seq.substr(it->start, it->length);
        //cout << s << endl;
    }
    //cout << "lengths: " << s.length() << ", " << p.length << endl;
    assert(s.length()==p.length);
    return s;
}

vector<LocalNode*> LocalPRG::nodes_along_path(const Path& p)
{
    vector<LocalNode*> v;
    // for each interval of the path
    for (deque<Interval>::const_iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
	//cout << "looking at interval " << *it << endl;
        // find the appropriate node of the prg
        for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
        {
            if ((it->end > n->second->pos.start and it->start < n->second->pos.end) or (it->start == n->second->pos.start and it->end == n->second->pos.end))
            {
		v.push_back(n->second);
		//cout << "found node " << *(n->second) << " so return vector size is now " << v.size() << endl;
            } else if (it->end < n->second->pos.start)
            {
		//cout << "gave up" << endl;
		break;
	    } // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
        }
    }
    return v;
}

vector<Interval> LocalPRG::splitBySite(const Interval& i)
{
    // Splits interval by next_site based on substring of seq in the interval
    //cout << "splitting by site " << next_site << " in interval " << i << endl;

    // Split first by var site
    vector<Interval> v;
    string::size_type k = i.start;
    string d = buff + to_string(next_site) + buff;
    string::size_type j = seq.find(d, k);
    while (j!=string::npos and j+d.size()<=i.end) {
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

    //cout << "after splitting by " << d << " got intervals " << v[0];
    assert(v[0].start >= i.start);
    for(uint32_t l=1; l!=v.size(); ++l)
    {
	//cout << " " << v[l];
        assert(v[l-1].end <= v[l].start|| assert_msg(v[l-1].end << ">" << v[l].start << " giving overlapping intervals  " << v[l-1] << " and " << v[l]));
    }
    //cout << endl;
    assert(v.back().end <= i.end);

    //assert(v.size()==1 or v.size()==3 || assert_msg(v.size() << "!=1,3 with first intervals " << v[0] << " and " << v[1]));

    // then split by var site + 1
    vector<Interval> w;
    d = buff + to_string(next_site+1) + buff;
    for(uint32_t l=0; l!=v.size(); ++l)
    {
	//cout << "splitting interval " << v[l] << " by " << d << endl;
	k = v[l].start;
	j = seq.find(d, k);
        while (j!=string::npos and j+d.size()<=v[l].end) {
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
    if (v.size() == w.size() && v.size() == 3)
    {
	cout << "There was something dodgy with var site " << next_site << ": found no separated alternates.";
	cout << " I'm going to assume for now that this is as a result of straggly ends of sequences which don't align nicely,";
	cout << " but you should check this. To handle, add an empty interval alternate." << endl;
	vector<Interval> x;
	for(uint32_t l=0; l!=w.size()-1; ++l)
	{
	    x.push_back(w[l]);
	}
	x.push_back(Interval(w[w.size()-2].end, w[w.size()-2].end));
	x.push_back(w[w.size()-1]);
        w = x;
    }
    //cout << "This was performed on sequence: " << seq << endl;
    //cout << "after splitting by " << d << " got intervals " << w[0];
    assert(w[0].start >= i.start);
    for(uint32_t l=1; l!=w.size(); ++l)
    {
	//cout << " " << w[l];
	assert(w[l-1].end <= w[l].start|| assert_msg(w[l-1].end << ">" << w[l].start << " giving overlapping intervals  " << w[l-1] << " and " << w[l] << " when splitting seq :" << seq.substr(i.start, i.length)));
    }
    //cout << endl;
    assert(w.back().end <= i.end);
    return w;
}

vector<uint32_t> LocalPRG::build_graph(const Interval& i, const vector<uint32_t>& from_ids, uint32_t current_level)
{
    // we will return the ids on the ends of any stretches of graph added
    vector<uint32_t> end_ids;

    // save the start id, so can add 0, and the last id to the index at level 0 at the end
    uint32_t start_id = next_id;

    // add nodes
    string s = seq.substr(i.start, i.length); //check length correct with this end...
    if (isalpha_string(s)) // should return true for empty string too
    {
	prg.add_node(next_id, s, i);
	//max_level = max(max_level, current_level);
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
	//cout << "split by site " << next_site << " in interval " << i << endl;
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
	//cout << "add pre site string " << s << endl;
	prg.add_node(next_id, s, v[0]);
	uint32_t pre_site_id = next_id;
	//max_level = max(max_level, current_level);
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
	    //cout << "add segment " << j << " /interval " << v[j] << endl;
	    vector<uint32_t> w = build_graph(v[j], mid_ids, current_level+1);
	    end_ids.insert(end_ids.end(), w.begin(), w.end());
	}
	// add this varsite to index
	prg.add_varsite (current_level + 1, pre_site_id, next_id);
	// add end interval
	end_ids = build_graph(v.back(), end_ids, current_level);
    }
    if (start_id == 0)
    {
	assert(end_ids.size()==1);
	prg.add_varsite (0, start_id, next_id - 1);
    }
    return end_ids;
}

void LocalPRG::minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k)
{
    //cout << "START SKETCH FUNCTION" << endl;
    vector<Path> walk_paths;
    deque<Interval> d;
    Path current_path, kmer_path;
    string kmer;
    pair<uint64_t, uint64_t> kh;
    uint64_t smallest = std::numeric_limits<uint64_t>::max();

    for (map<uint32_t,LocalNode*>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	//cout << "Processing node" << it->second->id << endl;
	
	// if part of the node has already been included in a mini, we start after this part
	for (uint32_t i=it->second->sketch_next; i!=max(it->second->sketch_next+w+k, it->second->pos.end+2)-w-k;)
        {
            //cout << "Node is long enough for whole w+k-1 length to fit in, " << it->second->pos.start << " <= " << i << " <= " << it->second->pos.end - w - k + 1 << endl; 
            assert(i <= it->second->pos.end - w - k + 1);
            assert(i >= it->second->pos.start);
	    
	    // note that as we are sketching disjoint kmers, have to work out from scratch for each window
            smallest = std::numeric_limits<uint64_t>::max();
            // find the lex smallest kmer in the window
            for (uint32_t j = 0; j < w; j++)
            {
                d = {Interval(i+j, i+j+k)};
                kmer_path.initialize(d);
                kmer = string_along_path(kmer_path);
                kh = kmerhash(kmer, k);
                smallest = min(smallest, min(kh.first, kh.second));
            }
            for (uint32_t j = 0; j < w; j++)
            {
                d = {Interval(i+j, i+j+k)};
                kmer_path.initialize(d);
                kmer = string_along_path(kmer_path);
                kh = kmerhash(kmer, k);
                if (kh.first == smallest)
                {
                    //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                    idx->add_record(kh.first, id, kmer_path, 0);
                    kmer_paths.push_back(kmer_path);
		    it->second->sketch_next = kmer_path.path.back().end;
                } else if (kh.second == smallest)
                {
                    //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                    idx->add_record(kh.second, id, kmer_path, 1);
                    kmer_paths.push_back(kmer_path);
		    it->second->sketch_next = kmer_path.path.back().end;
                }
            }
	    i=it->second->sketch_next;
	    if (i+w+k > it->second->pos.end+1)
	    {
		break;
	    } // too close to the end of this node for this method
	}

	// now sketch the remaining bit which overlaps the end boundary. If we have an empty node, still want a kmer covering, 
	// with corresponding empty interval at start. 
        for (uint32_t i=it->second->sketch_next; i<max(it->second->pos.end, it->second->pos.start+1);)
        {
	    //cout << "start " << it->second->pos.start << " <= " << i << " <= " << it->second->pos.end << " end" << endl;
	    if (i>=max(it->second->pos.end, it->second->pos.start+1))
	    {
		//cout << "I don't know how I got here" << endl;
		break;
	    }
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
                        smallest = min(smallest, min(kh.first, kh.second));
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
                        if (kh.first == smallest)
                        {
                            //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                            idx->add_record(kh.first, id, kmer_path, 0);
                            kmer_paths.push_back(kmer_path);
			    vector<LocalNode*> n = nodes_along_path(kmer_path);
			    for (vector<LocalNode*>::iterator it3 = n.begin(); it3!=n.end(); ++it3)
			    {
				//cout << "for node " << (*it3)->id << " sketch next was " << (*it3)->sketch_next << endl;
				//cout << " and will now be " << min((*it3)->pos.end, kmer_path.path.back().end) << endl;
				(*it3)->sketch_next = min((*it3)->pos.end, kmer_path.path.back().end);
			    }
                        } else if (kh.second == smallest)
                        {
                            //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                            idx->add_record(kh.second, id, kmer_path, 1);
                            kmer_paths.push_back(kmer_path);
			    vector<LocalNode*> n = nodes_along_path(kmer_path);
                            for (vector<LocalNode*>::iterator it3 = n.begin(); it3!=n.end(); ++it3)
                            {   
				//cout << "for node " << (*it3)->id << " sketch next was " << (*it3)->sketch_next << endl;
                                //cout << " and will now be " << min((*it3)->pos.end, kmer_path.path.back().end) << endl;
                                (*it3)->sketch_next = min((*it3)->pos.end, kmer_path.path.back().end);
                            }
                        }
		    }
		}
	    }
	    i=it->second->sketch_next;
	    //cout << "next i: " << i << " and it->second->pos.end: " << it->second->pos.end << endl;
            if (i+w+k > it->second->pos.end+1)
            {
		break;
	    } else { // too close to the end of this node for this method
		//cout << i << "<=" << it->second->pos.end-w-k+1 << endl;
	    }
	}
    }
    return;
}

/*void LocalPRG::minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k)
{
    //cout << "START SKETCH FUNCTION" << endl;
    vector<Path> walk_paths;
    deque<Interval> d;
    Path current_path, kmer_path, prev_path;
    string kmer;
    pair<uint64_t, uint64_t> kh;
    uint64_t smallest = std::numeric_limits<uint64_t>::max();

    for (map<uint32_t,LocalNode*>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	//cout << "Processing node" << it->second->id << endl;
	// if node has long seq, there will be no branching until  reach len-(w+k-1)th position
	for (uint32_t i=it->second->pos.start; i!=max(it->second->pos.start+w+k, it->second->pos.end+1)-w-k; ++i)
	{
	    //cout << "Node has been deemed long enough and " << it->second->pos.start << " <= " << i << " <= " << it->second->pos.end - w - k + 1 << endl; 
	    assert(i < it->second->pos.end - w - k + 1);
      	    assert(i >= it->second->pos.start);
            // if window pos is first for node or the previous minimizer starts outside current window, calculate new mini from scratch
            if((i == it->second->pos.start) or (prev_path.start < i))
            {
                smallest = std::numeric_limits<uint64_t>::max();
                // find the lex smallest kmer in the window
                for (uint32_t j = 0; j < w; j++)
                {
		    d = {Interval(i+j, i+j+k)};
		    kmer_path.initialize(d);
                    kmer = string_along_path(kmer_path);
                    kh = kmerhash(kmer, k);
                    smallest = min(smallest, min(kh.first, kh.second));
                }
                for (uint32_t j = 0; j < w; j++)
                {
		    d = {Interval(i+j, i+j+k)};
                    kmer_path.initialize(d);
                    kmer = string_along_path(kmer_path);
                    kh = kmerhash(kmer, k);
                    if (kh.first == smallest)
                    {
			//cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                        idx->add_record(kh.first, id, kmer_path, 0);
			kmer_paths.push_back(kmer_path);
			update_minimizer_counts_for_nodes(kmer_path);
                        prev_path = kmer_path;
                    } else if (kh.second == smallest)
                    {
                        //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                        idx->add_record(kh.second, id, kmer_path, 1);
			kmer_paths.push_back(kmer_path);
			update_minimizer_counts_for_nodes(kmer_path);
                        prev_path = kmer_path;
                    }
                }
            } else {
            // otherwise only need to do something if the kmer from the newest position at end of new window is smaller or equal to the previous smallest
		d = {Interval(i+w-1, i+w-1+k)};
                kmer_path.initialize(d);
                kmer = string_along_path(kmer_path);
                kh = kmerhash(kmer, k);
                //cout << "Last kh for wpos: " << kh << " compared to previous smallest: " << smallest << endl;
                if(kh.first <= smallest)
                {
		    //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                    idx->add_record(kh.first, id, kmer_path, 0);
		    kmer_paths.push_back(kmer_path);
		    update_minimizer_counts_for_nodes(kmer_path);
                    prev_path = kmer_path;
                } else if(kh.second <= smallest)
                {
                    //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                    idx->add_record(kh.second, id, kmer_path, 1);
		    kmer_paths.push_back(kmer_path);
		    update_minimizer_counts_for_nodes(kmer_path);
                    prev_path = kmer_path;
                }
                smallest = min(smallest, min(kh.first, kh.second));
            }
        }

	
	for (uint32_t i=max(it->second->pos.start+w+k, it->second->pos.end+1)-w-k; i!=it->second->pos.end; ++i)
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
		        smallest = min(smallest, min(kh.first, kh.second));
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
		    	if (kh.first == smallest)
		    	{
                            //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
			    idx->add_record(kh.first, id, kmer_path, 0);
			    kmer_paths.push_back(kmer_path);
			    update_minimizer_counts_for_nodes(kmer_path);
		    	} else if (kh.second == smallest)
                        {
                            //cout << "add record: " << kmer << " " << id << " " << kmer_path << endl;
                            idx->add_record(kh.second, id, kmer_path, 1);
			    kmer_paths.push_back(kmer_path);
			    update_minimizer_counts_for_nodes(kmer_path);
                        }
		    }
                }
	    }
	}
    }
}*/

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

/*void LocalPRG::get_covgs(MinimizerHits* minimizer_hits)
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
}*/

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

/*void LocalPRG::update_covg_with_hits(deque<MinimizerHit*>& mhs)
{
    // iterate over map first, because that is slow bit
    for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
    {
	// for each hit
        for (deque<MinimizerHit*>::iterator mh = mhs.begin(); mh != mhs.end(); ++mh)
	{
	    //for each interval of hit
	    for (deque<Interval>::const_iterator it=(*mh)->prg_path.path.begin(); it!=(*mh)->prg_path.path.end(); ++it)
	    {
                if (it->end > n->second->pos.start and it->start < n->second->pos.end)
                {
		    //cout << "added covg" << endl;
                    n->second->covg += min(it->end, n->second->pos.end) - max(it->start, n->second->pos.start);
                } //else if (it->end < n->second->pos.start)
                //{ break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
            }
	}
    }
}*/

/*void LocalPRG::update_minimizer_counts_for_nodes(Path& p)
{
    // update total num_minis for PRG
    //num_minis+=1;

    // then for each interval of the record...
    for (deque<Interval>::const_iterator it=p.path.begin(); it!=p.path.end(); ++it)
    {
        // ...update the num_minis counts on the appropriate node(s) of the prg
        for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
        {
            if (it->end > n->second->pos.start and it->start < n->second->pos.end)
            {
                n->second->num_minis += 1;
            } else if (it->end < n->second->pos.start)
            { break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
        }
    }
}*/
void LocalPRG::update_kmers_on_node_paths(vector<MaxPath>& vmp)
{
    for (uint direction=0; direction!=3; ++direction)
    {
	update_kmers_on_node_path(vmp[direction], kmer_path_probs[direction]);
    }
    return;
}

void LocalPRG::update_kmers_on_node_path(MaxPath& mp, const vector<float>& kp_probs)
{
    assert(mp.kmers_on_path.size() == kmer_paths.size());
    
    deque<Interval>::const_iterator it;
    vector<LocalNode*>::const_iterator node;
    bool found, reject, added_m_or_n, found_branch;
    vector<uint32_t> nums;
    set<uint32_t> subnums;

    cout << "looking for kmers overlapping node_path " << endl;

    for (uint32_t n=0; n!=mp.npath.size(); ++n)
    {
	cout << *(mp.npath[n]) << " ";
    }
    cout << endl; 

    // find the set of kmer_paths which overlap this node_path
    for (uint32_t n=0; n!=kmer_paths.size(); ++n)
    {
        it=kmer_paths[n].path.begin();
        node=mp.npath.begin();
	found = false;
	reject = false;
	
        if (mp.kmers_on_path[n] == true)
	{
	    found = true;
	}

        while (reject == false and found == false and node!=mp.npath.end() and it!=kmer_paths[n].path.end())
	{
            if ((it->end > (*node)->pos.start and it->start < (*node)->pos.end) or (*it == (*node)->pos))
            {
                // then this node overlaps this interval of the kmer
                // it needs to continue to overlap to end of kmer or node_path
                found = true;
		//cout << endl << "found start of a match: " <<  *it << " " << (*node)->pos << endl;
		//cout << "note " << (reject == false) << " " << (node!=mp.npath.end()) << " " << (it!=kmer_paths[n].path.end()) << endl; 
		while (reject == false and node!=mp.npath.end() and it!=kmer_paths[n].path.end())
		{
		    if ((it->end > (*node)->pos.start and it->start < (*node)->pos.end) or (*it == (*node)->pos))
                    {
			//cout << *it << " " << (*node)->pos << " match" << endl;
			node++;
			it++;
		    } else {
			// we have stopped matching and not because we reached the end of the kmer or node path
			//cout << "reject: " << kmer_paths[n] << " since " << *it << " and " << (*node)->pos  << " do not match " << endl;
			reject = true;
		    }
		}
		//cout << "end of while" << endl;
	    } else {
		// no match for this node and inteval
		// a match has to start either with the first node of node_path, or first interval of kmer
		// try iterating through combinations
		if (node==mp.npath.begin() and it!=kmer_paths[n].path.end())
		{
		    it++;
		} else {
		    it=kmer_paths[n].path.begin();
		    node++;
		}
	    }
	}

        // now if it was found and not rejected, add to nums;
        if (found == true and reject == false)
	{
            nums.push_back(n);
	    cout << "found kmer match for " << kmer_paths[n] << endl;
	}
    }

    cout << "have a list of overlapping kmers, now chose between ones which branch" << endl;

    // if have several covering directions at a branching/closing point, want only one of the options so take a subset
    for (uint32_t n=0; n!=nums.size(); ++n)
    {
        // if we haven't already decided to eliminate the kmer because it represents an alternative branch path, find kmers that branch from this kmer
        if (subnums.find(n) == subnums.end())
        {
            added_m_or_n = false;
            found_branch = false;
            for (uint32_t m=n+1; m!=nums.size(); ++m)
            {
                // for any pair of numbers, if they represent branching paths, add the one with the smaller probability into the subset
                //cout << n << " vs " << m << endl;
                if (kmer_paths[nums[n]].is_branching(kmer_paths[nums[m]]))
                {
                    found_branch = true;
                    if (kp_probs[nums[n]] > kp_probs[nums[m]])
                    {
                        cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch but prob " << kp_probs[nums[n]] << " > " << kp_probs[nums[m]] << endl;
                        subnums.insert(m);
                        added_m_or_n = true;    //we have kept/added n
                    } else if (kp_probs[nums[m]] > kp_probs[nums[n]]) {
                        cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch but prob " << kp_probs[nums[n]] << " < " << kp_probs[nums[m]] << endl;
                        subnums.insert(n);
                        added_m_or_n = true; //we have kept/added m
                    } else {
                        cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch and have same prob" << endl;
                        subnums.insert(m);
                        subnums.insert(n);
                    }
                } else {
                    cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " do not branch" << endl;
                }
            }
            // if we haven't kept any of the branching options, remove n from the subnums set
            if (found_branch == true and added_m_or_n == false)
            {
                cout << "found many branching kmers with same prob, so add " << (kmer_paths[nums[n]]) << endl;
                subnums.erase(n);
            }
        }
    }

    // now for all numbers in nums not in subnums, indicate to keep
    for (uint32_t n=0; n!=nums.size(); ++n)
    {
        if (std::find(subnums.begin(), subnums.end(), n) == subnums.end())
        {
            //cout << "Found a kmer " << kmer_paths[nums[n]] << " which branches from node " << node->id << " " << node->pos << " and which has prob " << kmer_path_probs[nums[n]] << endl;
            mp.kmers_on_path[nums[n]] = true;
        } else {
	    mp.kmers_on_path[nums[n]] = false;
	}
    }
    cout << "done identifying kmers on path" << endl;
    //assert(std::accumulate(mp.kmers_on_path.begin(), mp.kmers_on_path.end(), 0) > 0);
}

void LocalPRG::get_kmer_path_hit_counts(const PanNode* pnode)
{
    // note that within the foundHits for a pnode, hits which map to same read, prg and strand will be grouped together
    // we want to know for each prg minimizer, how much support there is from the reads
    // count how many hits against each minimizing_kmer of prg
    vector<uint32_t> counts(kmer_paths.size(), 0);
    kmer_path_hit_counts.resize(3, counts);
    cout << "there are " << kmer_paths.size() << " minimizing kmers to count hits against" << endl;
    set<MinimizerHit*, pComp_path>::iterator mh_previous = pnode->foundHits.begin();
    uint32_t num_hits = 1, num_kmers=0;
    for (set<MinimizerHit*, pComp_path>::iterator mh = (pnode->foundHits.begin()); mh != pnode->foundHits.end();)
    {
        mh++;
        if ((mh==pnode->foundHits.end()) or !((*mh)->prg_path == (*mh_previous)->prg_path) or ((*mh)->strand != (*mh_previous)->strand))
        {
            for (uint32_t i=num_kmers; i!=kmer_paths.size(); ++i)
            {
                if (kmer_paths[i]==(*mh_previous)->prg_path)
                {
                    //cout << "found path " << kmer_paths[i] << " == " << (*mh_previous)->prg_path << " with hit count " << num_hits << endl;
                    kmer_path_hit_counts[(*mh_previous)->strand][i] = num_hits;
		    kmer_path_hit_counts[2][i] += num_hits; // vector for both forward and reverse hits to catch any partially inverted copies
                    //cout << "kmer_path_hit_counts[" << i << "] = " << num_hits << endl;
                    num_hits = 1;
                    num_kmers +=1;
                    break;
                }
            }
            assert(num_hits==1);
        } else {
            num_hits += 1;
            //cout << num_hits << " ";
        }
        mh_previous++;
    }
    // check all the hits have been added by this process
    uint32_t sum_of_elems = std::accumulate(kmer_path_hit_counts[2].begin(), kmer_path_hit_counts[2].end(), 0);
    cout << "after adding counts of the " << pnode->foundHits.size() << " hits, have still got " << sum_of_elems << " hits added to tallies" << endl;
    assert(sum_of_elems==pnode->foundHits.size());
    return;
}

void LocalPRG::get_kmer_path_probs(const PanNode* pnode, uint32_t k, float e_rate)
{
    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    cout << "work out prob of seeing this number of hits against a kmer assuming it is truly present" << endl;
    float p_kmer, p_max, p_min, p=1/exp(e_rate*k);
    uint32_t n = pnode->foundReads.size(), big_p_count=0;
    cout << "n: " << n << ", p: " << p << endl;
    cout << "count:0 " << nchoosek(n,0) << " * " << pow(p,0) << " * " << pow(1-p,n-0) << endl;
    cout << "count:1 " << nchoosek(n,1) << " * " << pow(p,1) << " * " << pow(1-p,n-1) << endl;
    for (uint32_t direction=0; direction!=3; ++direction) // directions 0,1,2 correspond to forward hit, rev_complement hit and either/both
    {	
	kmer_path_probs.push_back({});
        p_max=numeric_limits<float>::lowest(), p_min=0;
        for (uint32_t i=0; i!=kmer_path_hit_counts[direction].size(); ++i)
        {
            p_kmer = log(nchoosek(n, kmer_path_hit_counts[direction][i])*pow(p,kmer_path_hit_counts[direction][i])*pow(1-p,n-kmer_path_hit_counts[direction][i]));
            cout << kmer_paths[i] << " " << kmer_path_hit_counts[direction][i] << " " << p_kmer << endl;
            kmer_path_probs[direction].push_back(p_kmer);
            p_max = max(p_max, p_kmer);
            p_min = min(p_min, p_kmer);
            if(p_kmer>-0.5)
            {
                //cout << p_kmer << " ";
                big_p_count+=1;
            }
        }
        cout << "for " << direction << " direction found " << big_p_count << " log probs > " << -0.5 << " with max and min log probs " << p_max << ", " << p_min << endl;
        cout << endl;
    }
    return;
}


void LocalPRG::infer_most_likely_prg_paths_for_corresponding_pannode(const PanNode* pnode, uint32_t k, float e_rate)
{
    // start by counting how many hits against each minimizing_kmer of prg
    cout << now() << "start by counting how many hits against each minimizing kmer in prg" << endl;
    get_kmer_path_hit_counts(pnode);

    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    cout << now() << "next work out prob of seeing this number of hits against a kmer assuming it is truly present" << endl;
    get_kmer_path_probs(pnode, k, e_rate);

    //now we iterate through the graph from the outmost level to the lowest level working out the most likely path(s)
    //need a data structure to remember what the most probable path(s) were for var sites already considered
    //the max_path_index, stored by the LocalPRG class

    vector<vector<MaxPath>> u, v, w; // w <- u <=> v
    vector<LocalNode*> x;
    vector<MaxPath> t; // each of u,v,w contains items of form t, with 3 components corresponding to fwd,rev,both
    vector<bool> y(kmer_paths.size(),false);
    float max_prob, max_mean_prob, next_largest;


    // start with the outmost level
    uint8_t max_level = 0;
    for (auto const& element : prg.index) {
        max_level = max(max_level, element.first);
    }

    // and for each level..
    for (uint level = max_level; level <= max_level; --level)
    {
	cout << endl << now() << "Looking at level " << level << endl;
        // ...for each varsite at this level...
        for (uint i = 0; i!=prg.index[level].size(); ++i)
        {
            // add the first node of each alternate allele for the varsite to a vector
            uint32_t pre_site_id = prg.index[level][i].first, post_site_id = prg.index[level][i].second;
	    if (level == 0)
	    {
		assert(pre_site_id == 0);
		cout << now() << "finally find best path through whole prg" << endl;
		x.push_back(prg.nodes[pre_site_id]);
		t.resize(3, MaxPath(x, y, 0));
		u.push_back(t);
		x.clear();
		t.clear();
	    } else {
	        cout << now() << "Looking at varsite number " << i << " at this level, from the outnodes of " << pre_site_id << " to the innodes of " << post_site_id << endl;
		// we want the index to be inclusive of first node/pre_site_id prob, but exclusive of end node prob
                for (uint j = 0; j!=prg.nodes[pre_site_id]->outNodes.size(); ++j)
                {
		    cout << now() << "add " << pre_site_id << "->" << prg.nodes[pre_site_id]->outNodes[j]->id << " to u" << endl;
		    x.push_back(prg.nodes[pre_site_id]);
                    x.push_back(prg.nodes[pre_site_id]->outNodes[j]);
		    cout << now() << "make maxpath" << endl;
		    t.resize(3, MaxPath(x, y, 0));
                    u.push_back(t);
                    x.clear();
                    t.clear();

		    update_kmers_on_node_paths(u.back());
                }
	    }

	    assert (u.size()>0); // we just added either start node, or branching outnodes of a varsite to it!

            // then until we reach the end varsite:
            while (u.size()>0)
            {
		cout << now() << "restart loop with u.size() == " << u.size() << endl;
                // for each tuple in u
                for (uint j = 0; j!=u.size(); ++j)
                {
		    assert(u[j].size() == 3 || assert_msg(u[j].size() << "== u[j].size()=/= 3")); // should have fwd, rev, both MaxPaths
		    assert(u[j][0].npath.back()->id == u[j][1].npath.back()->id); // the rest of the path may differ
		    assert(u[j][0].npath.back()->id == u[j][2].npath.back()->id);
                    // if the path in the tuple ends at the end of the varsite, it is a done path, add to w
                    if (u[j][0].npath.back()->id == post_site_id)
                    {
			cout << now() << "finished paths: " << endl;
			for (uint32_t dir = 0; dir != 3; ++dir)
			{
			    cout << "direction-" << dir << " ";
			    for (uint32_t m = 0; m!= u[j][dir].npath.size(); ++m)
			    {
			        cout << u[j][dir].npath[m]->id << "->";
			    }
			    cout << "p: " << u[j][dir].get_prob(kmer_path_probs[dir]) << endl;
			}
                        w.push_back(u[j]);
                    } else {
                        // otherwise look up the last node in the max_path_index
                        // if it is there, then we have a path to extend by, 
                        // and extend the path for each of the maximal paths through the sub_varsite
                        // and add the updated path/prob pairs to v
                        cout << now() << "extend paths: " << endl;
			for (uint32_t dir = 0; dir != 3; ++dir)
			{
			    cout << "direction-" << dir << " ";
			    for (uint32_t m = 0; m!= u[j][dir].npath.size(); ++m)
                            {
                                cout << u[j][dir].npath[m]->id << "->";
                            }
                            cout << "p: " << u[j][dir].get_prob(kmer_path_probs[dir]) << endl;
			}
                        map<uint32_t, vector<vector<MaxPath>>>::iterator it=max_path_index.find(u[j][0].npath.back()->id);
                        if (it != max_path_index.end())
                        {
                            //extend node path with the max paths seen before from this point
                            cout << now() << "found " << u[j][0].npath.back()->id << " in max_path_index" << endl;
			    assert(it->second.size() == 3);
			    for (uint32_t dir = 0; dir != 3; ++dir)
			    {
				assert(it->second[dir].size() == 1);
			        u[j][dir].extend(it->second[dir][0]);
			        cout << now() << "extended npath with found npath giving: ";
			        for (uint32_t m = 0; m!= u[j][dir].npath.size(); ++m)
                                {
                                    cout << u[j][dir].npath[m]->id << "->";
                                }
			        cout << endl;
			    }
			} else {
			    //otherwise extend with the outnode of the last node in node_path
			    cout << now() << "did not find " << u[j][0].npath.back()->id << " in max_path_index, so add outnode" << endl;
			    assert(u[j][0].npath.back()->outNodes.size() == 1); // if the back is the start of a bubble, should already be in index!
			    for (uint32_t dir = 0; dir != 3; ++dir)
			    {
			        u[j][dir].npath.push_back(u[j][dir].npath.back()->outNodes[0]);
			        u[j][dir].kmers_on_path = y; 
			    }
			}
			// either way, add extended node_path to v
			update_kmers_on_node_paths(u[j]);
			v.push_back(u[j]);
                    }
                }
                // once done for all of what was in u, set u = v
                u = v;
		v.clear();
            }

	    assert(w.size()>0); // we have to have found some paths
	    assert(u.size()==0);
	    assert(t.size()==0);

	    // at this point, all the paths should have up to date vectors of bools representing which kmers overlapped, 
	    // including the pre_site_id node and the post_site_id
	    // Can now work out the max of the mean probs along path for each direction

	    for (uint32_t dir = 0; dir != 3; ++dir)
	    {
                max_prob = numeric_limits<float>::lowest(); 
		max_mean_prob = numeric_limits<float>::lowest();
		next_largest = numeric_limits<float>::lowest();

                // find max for all probs we could work out values for
                cout << now() << "Finding max of: " << endl;
                for (uint n = 0; n!=w.size(); ++n)
                {
		    if (w[n][dir].get_mean_prob(kmer_path_probs[dir]) > max_mean_prob)
		    {
		        next_largest = max_mean_prob;
                        max_mean_prob = w[n][dir].mean_prob;
		    }
		    cout << w[n][dir].mean_prob << ", ";
                }
	        cout << endl;
	        // for ties, use max_prob
	        for (uint n = 0; n!=w.size(); ++n)
                {
                    if (w[n][dir].mean_prob == max_mean_prob)
                    {
                        max_prob = max(max_prob, w[n][dir].get_prob(kmer_path_probs[dir]));
                    }
                }
	        cout << now() << "max_prob (mean) for paths at this varsite: " << max_mean_prob << " and max_prob: " << max_prob << endl;

                // now add (a) path achieving max to max_path_index
	        // if there are multiple such paths, we just add the first
	        // note that in the case pre_site_id == 0 and level == 0, 
	        // we may overwrite a previous entry to the index
                for (uint n = 0; n!=w.size(); ++n)
                {
                    if ((max_mean_prob == numeric_limits<float>::lowest() and w[n][dir].mean_prob == 0) or (w[n][dir].mean_prob == max_mean_prob and w[n][dir].prob == max_prob))
                    {
			t.push_back(w[n][dir]);
		        cout << now() << "Add path to index: ";
		        for (uint32_t m = 0; m!= w[n][dir].npath.size(); ++m)
		        {
			    cout << w[n][dir].npath[m]->id << "->";
		        }
		        cout << "mean_p: " << w[n][dir].mean_prob << " compared to next largest " << next_largest << ", and p: " << w[n][dir].prob << endl;
		        break;
                    }	
                }
		assert(t.size() == 1);
		u.push_back(t);
		t.clear();
	    }

	    assert(u.size() == 3);
            max_path_index[pre_site_id] = u;
	    w.clear();
	    u.clear();
        }
    }

    // to finish, label the paths by which direction they came from, and sort so
    // most likely comes first
    max_path_index[0][0][0].direction = "forward";
    max_path_index[0][1][0].direction = "reverse-complement";
    max_path_index[0][2][0].direction = "both/inversion";
    sort( max_path_index[0].begin(), max_path_index[0].end(), VMPgreater() );
    return;
}

void LocalPRG::write_max_paths_to_fasta(const string& filepath)
{
    assert(max_path_index.size()>0);
    map<uint32_t, vector<vector<MaxPath>>>::iterator it=max_path_index.find(0);
    assert(it!=max_path_index.end());
    assert(max_path_index[0].size()==3);
    assert(max_path_index[0][0].size()>0);
    assert(max_path_index[0][1].size()>0);
    assert(max_path_index[0][2].size()>0);

    ofstream handle;
    handle.open (filepath);
    for (uint i = 0; i!= 3; ++i)
    {
        handle << ">" << name << "." << max_path_index[0][i][0].direction << "\t P(data|sequence)=" << max_path_index[0][i][0].prob << endl;
	for (uint j = 0; j!= max_path_index[0][i][0].npath.size(); ++j)
	{
            handle << max_path_index[0][i][0].npath[j]->seq;
        }
        handle << endl;

	/*// and the reverse complement
	handle << ">rc_" << id << "." << i << "\t P(data|sequence)=" << max_path_index[0][i].prob << endl;
        for (uint j = max_path_index[0][i][0].npath.size(); j!= 0; --j)
        {
            handle << rev_complement(max_path_index[0][i].npath[j-1]->seq);
        }
        handle << endl;*/
    }
    handle.close();
    return;
}
