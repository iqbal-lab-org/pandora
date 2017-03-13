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

using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p): next_id(0), buff(" "), next_site(5), id(i), name(n), seq(p)
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

    kmer_paths.reserve(50000);
    kmer_path_clashes.reserve(50000);
}

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

/*vector<LocalNode*> LocalPRG::nodes_along_string(const string& query_string)
{
    vector<vector<LocalNode*>> u,v;   // u <=> v
				      // ie reject paths in u, or extend and add to v
				      // then set u=v and continue
    u.reserve(100);
    v.reserve(100);
    vector<LocalNode*> npath;
    string candidate_string = "";
    
    assert(prg.nodes.size()>0); //otherwise empty prg -> segfault
    u = {{prg.nodes[0]}};

    while (u.size() > 0)
    {
	for (uint32_t i=0; i!=u.size(); ++i)
	{
	    for (uint32_t j=0; j!=u[i].size(); ++j)
	    {
		candidate_string += u[i][j]->seq;
	    }

	    for (uint32_t j=0; j!=u[i].back()->outNodes.size(); ++j)
	    {
		// if the start of query_string matches extended candidate_string, want to query candidate path extensions
		if ( query_string.substr(0,candidate_string.size()+u[i].back()->outNodes[j]->seq.size()) == candidate_string+u[i].back()->outNodes[j]->seq)
		{
		    if (candidate_string.size()+u[i].back()->outNodes[j]->seq.size() >= query_string.size())
		    {
			// we have now found the whole of the query_string
			u[i].push_back(u[i].back()->outNodes[j]);
			return u[i];
		    } else {
		        v.push_back(u[i]);
		        v.back().push_back(u[i].back()->outNodes[j]);
		    }
		}
	    }
    	    candidate_string = "";
	}
	u = v;
	v.clear();
    }
    // found no successful path, so return an empty vector
    return npath;
}*/

vector<LocalNode*> LocalPRG::nodes_along_path(const Path& p)
{
    vector<LocalNode*> v;
    v.reserve(100);
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

vector<Interval> LocalPRG::split_by_site(const Interval& i)
{
    // Splits interval by next_site based on substring of seq in the interval
    //cout << "splitting by site " << next_site << " in interval " << i << endl;

    // Split first by var site
    vector<Interval> v;
    v.reserve(4);
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
    w.reserve(20);
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
    end_ids.reserve(20);

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
	vector<Interval> v = split_by_site(i); // should have length at least 4
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
        mid_ids.reserve(20);
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
    cout << now() << "Sketch PRG " << name << " which has " << prg.nodes.size() << " nodes" << endl;

    // clean up after any previous runs
    kmer_prg.clear();
    kmer_paths.clear();
    for (map<uint32_t,LocalNode*>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	it->second->prev_kmer_paths.clear();
	it->second->sketch_next = it->second->pos.start;
    }

    // declare variables
    vector<Path> walk_paths;
    set<Path> ends;
    vector<LocalNode*> n;
    deque<Interval> d;
    Path current_path, kmer_path;
    string kmer;
    pair<uint64_t, uint64_t> kh;
    KmerHash hash;
    uint64_t smallest = std::numeric_limits<uint64_t>::max();

    // create a null start node in the kmer graph
    d = {Interval(0,0)};
    kmer_path.initialize(d);
    prg.nodes.begin()->second->prev_kmer_paths.insert(kmer_path);
    kmer_prg.add_node(kmer_path);

    // force every start kmer to be in the graph
    walk_paths = prg.walk(prg.nodes.begin()->second->id, 0, k);
    //assert(walk_paths.size()>=1 or prg.nodes.size()<=1);
    if (walk_paths.size() == 1 and walk_paths.back().path.size() == 1)
    {
        prg.nodes.begin()->second->prev_kmer_paths.clear();
    }
    for (uint i=0; i!=walk_paths.size(); ++i)
    {
        kmer = string_along_path(walk_paths[i]);
        kh = hash.kmerhash(kmer, k);

	// add to index, kmer_prg and kmer_paths
	idx->add_record(min(kh.first, kh.second), id, walk_paths[i], (kh.first>kh.second));
	kmer_prg.add_node(walk_paths[i]);
        kmer_prg.add_edge(kmer_prg.nodes[0]->path, walk_paths[i]);
	n = nodes_along_path(walk_paths[i]);
	for (uint m=0; m!=n.size(); ++m)
        {
            n[m]->sketch_next = min(n[m]->pos.end, n[m]->pos.start+1);
            n[m]->prev_kmer_paths.insert(walk_paths[i]);
	    n[m]->skip = true; 
        }
        if (walk_paths[i].end == n.back()->pos.end and n.back()->outNodes.size() > 0)
        {
            for (uint m=0; m!=n.back()->outNodes.size(); ++m)
            {
                n.back()->outNodes[m]->prev_kmer_paths.insert(walk_paths[i]);
            }
        }
        kmer_paths.push_back(walk_paths[i]);
    }

    // now process each of the LocalNodes in turn
    for (map<uint32_t,LocalNode*>::iterator it=prg.nodes.begin(); it!=prg.nodes.end(); ++it)
    {
	// if we have a null node which has already been sketched, skip
        while (it!=prg.nodes.end() and it->second->skip == true and it->second->pos.length == 0)
        {
            cout << "skip node " << *(it->second) << endl;
	    it++;
        }
	if (it==prg.nodes.end())
        {
            break;
        }

        cout << "Processing node" << *(it->second) << endl;
        smallest = std::numeric_limits<uint64_t>::max();

	// if part of the node has already been included in a mini, we start after this part
	for (uint32_t i=it->second->sketch_next; i<max(it->second->pos.end, it->second->pos.start+1);)
        {
	    //cout << "i: " << i << endl;
	    assert(i<max(it->second->pos.end, it->second->pos.start+1));

	    // if this is not the first minikmer for this node, and doesn't overlap
	    // a branch point, then only need to consider the newest kmer as the window
	    // shifts a base, and compare to old smallest, provided window still includes the old one
	    if (smallest != std::numeric_limits<uint64_t>::max() and i+w+k-1<it->second->pos.end and it->second->prev_kmer_paths.begin()->start >= i)
	    {
	        d = {Interval(i+w-1, i+w-1+k)};
                kmer_path.initialize(d);
                kmer = string_along_path(kmer_path);
                kh = hash.kmerhash(kmer, k);
		//cout << "compare last kmer " << min(kh.first, kh.second) << " with smallest kmer " << smallest << endl;
                if(kh.first <= smallest or kh.second <= smallest)
                {
		    idx->add_record(min(kh.first, kh.second), id, kmer_path, (kh.first>kh.second));
		    kmer_paths.push_back(kmer_path);
                    kmer_prg.add_node(kmer_path);
		    assert(it->second->prev_kmer_paths.size()==1);
		    kmer_prg.add_edge(*(it->second->prev_kmer_paths.begin()), kmer_path);
                    it->second->prev_kmer_paths.clear();
		    it->second->prev_kmer_paths.insert(kmer_path);
		    smallest = min(kh.first, kh.second);
		}
	    } else {

		//cout << "walk and minimize window" << endl;
                walk_paths = prg.walk(it->second->id, i, w+k-1);
	        //assert(walk_paths.size()>=1);
	        //cout << "found " << walk_paths.size() << " walk paths" << endl;
	        if (walk_paths.size() == 0)
	        {
		    break;
	        }

                for (vector<Path>::iterator it2=walk_paths.begin(); it2!=walk_paths.end(); ++it2)
                {
		    //cout << "minimize walk path " << *it2 << endl;
                    // find minimizer for this path 
                    smallest = std::numeric_limits<uint64_t>::max();
                    for (uint32_t j = 0; j != w; j++)
                    {
                        kmer_path = it2->subpath(i+j,k);
                        if (kmer_path.path.size() > 0)
                        {
                            kmer = string_along_path(kmer_path);
                            kh = hash.kmerhash(kmer, k);
                            smallest = min(smallest, min(kh.first, kh.second));
                        }
                    }
                    for (uint32_t j = 0; j != w; j++)
                    {
                        kmer_path = it2->subpath(i+j,k);
                        if (kmer_path.path.size() > 0)
                        {
                            kmer = string_along_path(kmer_path);
                            kh = hash.kmerhash(kmer, k);

			    n = nodes_along_path(kmer_path);
                            // add null nodes to end end of the kmer path if there is not enough length
                            // left to end of PRG for them to be included in another kmer
                            if (prg.walk(n.back()->id, n.back()->pos.end, w+k-1).size() == 0)
                            {
                                while (kmer_path.end >= n.back()->pos.end and n.back()->outNodes.size() == 1 and n.back()->outNodes[0]->pos.length == 0)
                                {
                                    kmer_path.add_end_interval(n.back()->outNodes[0]->pos);
                                    n.push_back(n.back()->outNodes[0]);
                                }
                            }
	
                            if ((kh.first == smallest or kh.second == smallest or n.back() == (--(prg.nodes.end()))->second) and 
				(find(kmer_paths.begin(), kmer_paths.end(), kmer_path)==kmer_paths.end()))
                            {
			        // add to index, kmer_prg and kmer_paths
			        idx->add_record(min(kh.first, kh.second), id, kmer_path, (kh.first>kh.second));
                                kmer_paths.push_back(kmer_path);
			        kmer_prg.add_node(kmer_path);
				assert(it->second->prev_kmer_paths.size() > 0);
				for (set<Path>::iterator l = it->second->prev_kmer_paths.begin(); l!=it->second->prev_kmer_paths.end(); ++l)
                                {
                                    if ((*l < kmer_path) and !(kmer_path.is_branching(*l)))
                                    {
                                        kmer_prg.add_edge(*l, kmer_path);
                                    }
                                }	
				assert(kmer_prg.nodes.back()->inNodes.size() > 0);
			        ends.insert(kmer_path);    
			        for (uint m=(it->second == n[0]); m<n.size(); ++m)
			        {
				    n[m]->prev_kmer_paths.insert(kmer_path);
				    n[m]->skip = true; // nb only used by null nodes, and only last kmers in prg can end on null kmers
			        }
				if (n.back()->pos.end == kmer_path.end)
				{
				    for (uint m=0; m!=n.back()->outNodes.size(); ++m)
                                    {
                                        n.back()->outNodes[m]->prev_kmer_paths.insert(kmer_path);
				    }
				}
			    } else if (kh.first == smallest or kh.second == smallest or n.back() == (--(prg.nodes.end()))->second) {
				ends.insert(kmer_path);
			    }
		        }
		    }
	        }
	        it->second->prev_kmer_paths = ends;
	        ends.clear();
	    }
	    i++;
            //cout << "sketch " << i << " next" << endl;
            if (i >= it->second->pos.end)
            {
                //cout << "finish node, as " << i << ">=" << it->second->pos.end << endl;
                break;
            }
	}
    }

    // now force adding of the last kmer for each path
    //cout << "add end nodes" << endl;
    walk_paths = prg.walk_back((--(prg.nodes.end()))->second->id, (--(prg.nodes.end()))->second->pos.end, k);
    for (uint i=0; i!=walk_paths.size(); ++i)
    {
	if (find((--(prg.nodes.end()))->second->prev_kmer_paths.begin(), (--(prg.nodes.end()))->second->prev_kmer_paths.end(), walk_paths[i]) == (--(prg.nodes.end()))->second->prev_kmer_paths.end())
	{
	    cout << "add previously unfound end kmer" << endl;
            kmer = string_along_path(walk_paths[i]);
            kh = hash.kmerhash(kmer, k);

            // add to index, kmer_prg and kmer_paths
            idx->add_record(min(kh.first, kh.second), id, walk_paths[i], (kh.first>kh.second));
	    n = nodes_along_path(walk_paths[i]);
	    kmer_prg.add_node(walk_paths[i]);
	    set<Path> s = n[0]->prev_kmer_paths;
	    cout << "prev_kmer_paths from n[0] " << *n[0] << " has size " << s.size() << endl;
	    set<Path> t;
	    cout << "while back " << *(kmer_prg.nodes.back()) << " has no innodes" << endl;
	    while (kmer_prg.nodes.back()->inNodes.size() == 0)
	    {
		cout << "try to find an innode from a set of size " << s.size() << endl;
		for (set<Path>::iterator l = s.begin(); l!=s.end(); ++l)
		{
		    cout << "consider adding edge to path " << *l << endl;
		    if ((*l < walk_paths[i]) and !(walk_paths[i].is_branching(*l)))
                    {
                        kmer_prg.add_edge(*l, walk_paths[i]);
		    } else if (walk_paths[i] < *l) {
			cout << "consider innodes instead" << endl;
		        t.insert(kmer_prg.get_innodes(*l).begin(), kmer_prg.get_innodes(*l).end());
		    } else {
			cout << walk_paths[i] << " and " << *l << " are not compatible" << endl;
		    }
		}
		s = t;
                t.clear();
	    }
	    assert(kmer_prg.nodes.back()->inNodes.size() > 0);
	    if (walk_paths.size()==0 and walk_paths[i].path.size() == 1)
	    {
	        n.back()->prev_kmer_paths.clear();
	    }
            kmer_paths.push_back(walk_paths[i]);
	}
	ends.insert(walk_paths[i]);
    }

    // create a null end node in the kmer graph
    //cout << "add null end" << endl;
    map<uint32_t,LocalNode*>::iterator it=--(prg.nodes.end());
    d = {Interval(it->second->pos.end,it->second->pos.end)};
    kmer_path.initialize(d);
    if (kmer_path == kmer_prg.nodes[0]->path)
    {
	return;
    }
    kmer_prg.add_node(kmer_path);

    for (set<Path>::iterator l = ends.begin(); l!=ends.end(); ++l)
    {
        kmer_prg.add_edge(*l, kmer_path);
    }

    // print, check and return
    //cout << "kmer prg: " << endl << kmer_prg << endl;
    kmer_prg.check(kmer_paths.size());
    return;
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

vector<LocalNode*> LocalPRG::localnode_path_from_kmernode_path(vector<KmerNode*> kmernode_path)
{
    cout << now() << "convert kmernode path to localnode path" << endl;
    vector<LocalNode*> localnode_path, kmernode;
    for (uint i=0; i!=kmernode_path.size(); ++i)
    {
	//cout << kmernode_path[i]->path << endl;
        kmernode = nodes_along_path(kmernode_path[i]->path);
	// if the start of the new localnode path is after the end of the previous, join up WLOG with top path
        while (localnode_path.size() > 0 and localnode_path.back()->outNodes.size() > 0 and kmernode[0]->id > localnode_path.back()->outNodes[0]->id)
        {
            localnode_path.push_back(localnode_path.back()->outNodes[0]);
        }
	// if new the localnodes in kmernode overlap old ones, truncate localnode_path
	while (localnode_path.size() > 0 and kmernode[0]->id <= localnode_path.back()->id)
	{
	    localnode_path.pop_back();
	}
	localnode_path.insert(localnode_path.end(), kmernode.begin(), kmernode.end());
	/*cout << "localnodepath: ";
	for (uint j=0; j!=localnode_path.size(); ++j)
	{
	    cout << *localnode_path[j] << "->";
	}*/
    }
    return localnode_path;
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
	    } else if (*it == n->second->pos) {
		n->second->covg += 1;
            } else if (it->end < n->second->pos.start)
            { break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
        }
    }

    // update the covg in the kmer_prg
    for (uint i=0; i!=kmer_prg.nodes.size(); ++i)
    {
	if (kmer_prg.nodes[i]->path == mh->prg_path)
	{
	    kmer_prg.nodes[i]->covg[mh->strand] += 1;
	    break;
	}
    }
}

void LocalPRG::update_kmers_on_node_paths(vector<MaxPath>& vmp)
{
    for (uint direction=0; direction!=3; ++direction)
    {
	update_kmers_on_node_path(vmp[direction], kmer_path_probs[direction]);
    }
    return;
}

vector<uint32_t> LocalPRG::find_overlapping_kmer_paths(MaxPath& mp)
{

    deque<Interval>::const_iterator it;
    vector<LocalNode*>::const_iterator node;

    bool found, reject;
    vector<uint32_t> nums;
    nums.reserve(40);

    // find the set of kmer_paths which overlap this node_path
    for (uint32_t n=0; n!=kmer_paths.size(); ++n)
    {
        if (mp.kmers_on_path[n] == 1)
        {
            nums.push_back(n);
        } else if (kmer_paths[n].end >= mp.npath[0]->pos.start and kmer_paths[n].start <= mp.npath.back()->pos.end) {
	    
            found = false;
            reject = false;

            it=kmer_paths[n].path.begin();
            node=mp.npath.begin();

            while (reject == false and found == false and node!=mp.npath.end())
            {
                if (it!=kmer_paths[n].path.end() and
                   ((it->end > (*node)->pos.start and it->start < (*node)->pos.end) or
                   (*it == (*node)->pos)))
                {
                    // then this node overlaps this interval of the kmer
                    // it needs to continue to overlap to end of kmer or node_path
                    found = true;
                    node++;
                    it++;
                    while (reject == false and node!=mp.npath.end() and it!=kmer_paths[n].path.end())
                    {
                        if (it!=kmer_paths[n].path.end() and
                            ((*it == (*node)->pos) or
			    (it->end > (*node)->pos.start and it->start < (*node)->pos.end)))
                        {
                            node++;
                            it++;
                        } else {
                            // we have stopped matching and not because we reached the end of the kmer or node path
                            reject = true;
			    break;
                        }
                    }
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
                mp.kmers_on_path[n] = 1;
                //cout << "found kmer match for " << kmer_paths[n] << endl;
            }
        }
    }
    return nums;
}

void LocalPRG::filter_branching_kmer_paths(MaxPath& mp, const vector<float>& kp_probs, const vector<uint32_t>& nums)
{
    bool added_m_or_n, found_branch;

    // if have several covering directions at a branching/closing point, want only one of the options so take a subset
    for (uint32_t n=0; n!=nums.size(); ++n)
    {
	added_m_or_n = false;
        found_branch = false;
        for (uint32_t m=n+1; m!=nums.size(); ++m)
	{
	    //check that these two kmers are not branching, removing one or both if they are
	    if ((mp.kmers_on_path[nums[n]] == 1 or found_branch == true) and mp.kmers_on_path[nums[m]] == 1)
	    {
                if (kmer_paths[nums[n]].is_branching(kmer_paths[nums[m]]))
                {
                    found_branch = true;
                    if (kp_probs[nums[n]] > kp_probs[nums[m]])
                    {
                        //cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch but prob " << kp_probs[nums[n]] << " > " << kp_probs[nums[m]] << endl;
                        mp.kmers_on_path[nums[m]] = 0;
                        added_m_or_n = true;    //we have kept/added n
                    } else if (kp_probs[nums[m]] > kp_probs[nums[n]]) {
                        //cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch but prob " << kp_probs[nums[n]] << " < " << kp_probs[nums[m]] << endl;
                        mp.kmers_on_path[nums[n]] = 0;
                        added_m_or_n = true; //we have kept/added m
                    } else {
                        //cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " branch and have same prob" << endl;
                        mp.kmers_on_path[nums[m]] = 0;
                        mp.kmers_on_path[nums[n]] = 0;
                    }
                } else {
		    //cout << "Kmers " << kmer_paths[nums[n]] << " and " << kmer_paths[nums[m]] << " do not branch" << endl;
	        }
	    }
	}
        // if we haven't kept any of the branching options, remove n from the subnums set
        if (found_branch == true and added_m_or_n == false)
        {
            //cout << "found many branching kmers with same prob, so add " << (kmer_paths[n]) << endl;
            mp.kmers_on_path[nums[n]] = 1;
        }
    }

    //cout << "done identifying kmers on path - found " << accumulate(mp.kmers_on_path.begin(), mp.kmers_on_path.end(), 0) << endl;
    //assert(std::accumulate(mp.kmers_on_path.begin(), mp.kmers_on_path.end(), 0) > 0);
}

void LocalPRG::update_kmers_on_node_path(MaxPath& mp, const vector<float>& kp_probs)
{
    assert(mp.kmers_on_path.size() == kmer_paths.size());
    assert(mp.kmers_on_path.size() == kp_probs.size());

    // find the set of kmer_paths which overlap this node_path
    vector<uint32_t> nums = find_overlapping_kmer_paths(mp);
    //cout << "have a list of " << nums.size() << " overlapping kmers, now chose between ones which branch" << endl;

    // filter branching kmer_paths
    filter_branching_kmer_paths(mp, kp_probs, nums);
    return;
}

void LocalPRG::get_kmer_path_hit_counts(const PanNode* pnode)
{
    // note that within the foundHits for a pnode, hits which map to same read, prg and strand will be grouped together
    // we want to know for each prg minimizer, how much support there is from the reads
    // count how many hits against each minimizing_kmer of prg
    vector<uint32_t> counts(kmer_paths.size(), 0);
    kmer_path_hit_counts.resize(3, counts);
    cout << now() << "There are " << kmer_paths.size() << " minimizing kmers to count hits against" << endl;
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
                    //num_kmers +=1;
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
    cout << now() << "After adding counts of the " << pnode->foundHits.size() << " hits, have got " << sum_of_elems << " hits added to tallies" << endl;
    assert(sum_of_elems==pnode->foundHits.size());
    return;
}

/*void LocalPRG::get_kmer_path_clash_counts()
{
    for (uint32_t i=0; i!=kmer_paths.size(); ++i)
    {
        kh = kmerhash(string_along_path(kmer_paths[i]), k);
	 
    vector<uint32_t> kmer_path_clashes; // number of places elsewhere in PRG where a kmer occurs TO DO!!

}*/

void LocalPRG::get_kmer_path_probs(const PanNode* pnode, uint32_t k, float e_rate)
{
    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    cout << now() << "For each PRG kmer, find log probability of seeing this number of hits assuming it is truly present" << endl;
    vector<float> probs(kmer_paths.size(), 0);
    kmer_path_probs.resize(3, probs);
    assert(kmer_path_probs.size() == 3);
    assert(kmer_path_probs[0].size() == kmer_paths.size());
    float p_kmer, p_max, p_min, p=1/exp(e_rate*k);
    uint32_t n = pnode->foundReads.size();
    cout << now() << "binomial parameters n: " << n << ", p: " << p << endl;
    for (uint32_t direction=0; direction!=3; ++direction) // directions 0,1,2 correspond to forward hit, rev_complement hit and either/both
    {	
        p_max=numeric_limits<float>::lowest(), p_min=0;
	assert(kmer_path_hit_counts[direction].size() == kmer_path_probs[direction].size());
        for (uint32_t i=0; i!=kmer_path_hit_counts[direction].size(); ++i)
        {
            p_kmer = lognchoosek(n, kmer_path_hit_counts[direction][i]) + kmer_path_hit_counts[direction][i]*log(p) + (n-kmer_path_hit_counts[direction][i])*log(1-p);
            kmer_path_probs[direction][i] = p_kmer;
            p_max = max(p_max, p_kmer);
            p_min = min(p_min, p_kmer);
        }
        cout << now() << "For " << direction << " direction found max and min log probs " << p_max << ", " << p_min << endl;
    }
    return;
}


void LocalPRG::infer_most_likely_prg_paths_for_corresponding_pannode(const PanNode* pnode, uint32_t k, float e_rate)
{
    // start by counting how many hits against each minimizing_kmer of prg
    get_kmer_path_hit_counts(pnode);

    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    get_kmer_path_probs(pnode, k, e_rate);

    //now we iterate through the graph from the outmost level to the lowest level working out the most likely path(s)
    //need a data structure to remember what the most probable path(s) were for var sites already considered
    //the max_path_index, stored by the LocalPRG class

    vector<vector<MaxPath>> u, v, w; // w <- u <=> v
    u.reserve(10);
    v.reserve(10);
    w.reserve(10);
    vector<LocalNode*> x;
    x.reserve(100);
    vector<MaxPath> t; // each of u,v,w contains items of form t, with 3 components corresponding to fwd,rev,both
    vector<int> y(kmer_paths.size(),0);
    float max_mean_prob;


    // start with the outmost level
    uint8_t max_level = prg.index.size() - 1;

    // and for each level..
    for (uint level = max_level; level <= max_level; --level)
    {
	cout << now() << "Looking at level " << level << endl;
        // ...for each varsite at this level...
        for (uint i = 0; i!=prg.index[level].size(); ++i)
        {
            // add the first node of each alternate allele for the varsite to a vector
            uint32_t pre_site_id = prg.index[level][i].first, post_site_id = prg.index[level][i].second;
	    assert(u.size()==0);
	    if (level == 0)
	    {
		assert(pre_site_id == 0);
		t.resize(3, MaxPath({prg.nodes[pre_site_id]}, y, 0));
		u = {t};
		t.clear();
	    } else {
		// we want the index to be inclusive of first node/pre_site_id prob, but exclusive of end node prob
                for (uint j = 0; j!=prg.nodes[pre_site_id]->outNodes.size(); ++j)
                {
		    t.resize(3, MaxPath({prg.nodes[pre_site_id], prg.nodes[pre_site_id]->outNodes[j]}, y, 0));
                    u.push_back(t);
		    t.clear();
		    update_kmers_on_node_paths(u.back());
                }
	    }
	
	    assert (u.size()>0); // we just added either start node, or branching outnodes of a varsite to it!

            // then until we reach the end varsite:
            while (u.size()>0)
            {
                // for each tuple in u
                for (uint j = 0; j!=u.size(); ++j)
                {
		    assert(u[j].size() == 3 || assert_msg(u[j].size() << "== u[j].size()=/= 3")); // should have fwd, rev, both MaxPaths
		    assert(u[j][0].npath.back()->id == u[j][1].npath.back()->id); // the rest of the path may differ
		    assert(u[j][0].npath.back()->id == u[j][2].npath.back()->id);
                    // if the path in the tuple ends at the end of the varsite, it is a done path, add to w
                    if (u[j][0].npath.back()->id == post_site_id)
                    {
                        w.push_back(u[j]);
                    } else {
                        // otherwise look up the last node in the max_path_index
                        // if it is there, then we have a path to extend by, 
                        // and extend the path for each of the maximal paths through the sub_varsite
                        // and add the updated path/prob pairs to v
                        map<uint32_t, vector<vector<MaxPath>>>::iterator it=max_path_index.find(u[j][0].npath.back()->id);
                        if (it != max_path_index.end())
                        {
                            //extend node path with the max paths seen before from this point
			    assert(it->second.size() == 3);
			    for (uint32_t dir = 0; dir != 3; ++dir)
			    {
				assert(it->second[dir].size() == 1);
			        u[j][dir].extend(it->second[dir][0]);
			    }
			} else {
			    //otherwise extend with the outnode of the last node in node_path
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
		max_mean_prob = numeric_limits<float>::lowest();
		int max_num_minis = 0;

                // find max for all probs we could work out values for
                for (uint n = 0; n!=w.size(); ++n)
                {
		    if (w[n][dir].get_mean_prob(kmer_path_probs[dir]) > max_mean_prob and w[n][dir].mean_prob != 0)
		    {
                        max_mean_prob = w[n][dir].mean_prob;
		    }
                }

		// for ties, pick the path with the most minimizers
		for (uint n = 0; n!=w.size(); ++n)
                {
                    if (w[n][dir].mean_prob == max_mean_prob)
                    {
			max_num_minis = max(max_num_minis, accumulate(w[n][dir].kmers_on_path.begin(), w[n][dir].kmers_on_path.end(), 0));
		    }
		}

                // now add (a) path achieving max to max_path_index
	        // if there are multiple such paths, we just add the first
	        // note that in the case pre_site_id == 0 and level == 0, 
	        // we may overwrite a previous entry to the index
	        assert(t.size() == 0);
                for (uint n = 0; n!=w.size(); ++n)
                {
                    if ((max_mean_prob == numeric_limits<float>::lowest() and w[n][dir].mean_prob == 0) or 
			(w[n][dir].mean_prob == max_mean_prob and accumulate(w[n][dir].kmers_on_path.begin(), w[n][dir].kmers_on_path.end(), 0) == max_num_minis))
                    {
			w[n][dir].get_prob(kmer_path_probs[dir]);
			t.push_back(w[n][dir]);
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

vector<float> LocalPRG::get_covered_maxpath_log_probs(uint dir, uint num_minis)
{
    //assert that we already have counts and probs
    assert(accumulate(kmer_path_hit_counts[2].begin(), kmer_path_hit_counts[2].end(), 0) > 0);

    //we iterate through the graph from the outmost level to the lowest level working out the covered path(s)
    vector<MaxPath> u, v, w; // w <- u <=> v
    map<uint32_t, vector<MaxPath>> bubble_paths;
				     // this time each has size 3, with an unknown number of options within that
    u.reserve(600000);
    v.reserve(600000);
    w.reserve(100000);
    vector<int> y(kmer_paths.size(),0);

    // start with the outmost level
    uint8_t max_level = prg.index.size() - 1;

    // and for each level..
    for (uint level = max_level; level <= max_level; --level)
    {
	cout << now() << "Looking at level " << level << endl;
        // ...for each varsite at this level...
        for (uint i = 0; i!=prg.index[level].size(); ++i)
        {
            // add the first node of each alternate allele for the varsite to a vector
            uint32_t pre_site_id = prg.index[level][i].first, post_site_id = prg.index[level][i].second;
	    assert(u.size()==0);
	    if (level == 0)
	    {
		assert(pre_site_id == 0);
		u.push_back(MaxPath({prg.nodes[pre_site_id]}, y, 0));
	    } else {
		// we want the index to be inclusive of first node/pre_site_id prob, but exclusive of end node prob
                for (uint j = 0; j!=prg.nodes[pre_site_id]->outNodes.size(); ++j)
                {
                    u.push_back(MaxPath({prg.nodes[pre_site_id], prg.nodes[pre_site_id]->outNodes[j]}, y, 0));
		    update_kmers_on_node_path(u.back(), kmer_path_probs[dir]);
                }
	    }
	
	    assert (u.size()>0); // we just added either start node, or branching outnodes of a varsite to it!

            // then until we reach the end varsite:
            while (u.size()>0)
            {
		cout << "u.size()==" << u.size() << endl;
                // for each path in u
                for (uint j = 0; j!=u.size(); ++j)
                {
                    // if the path ends at the end of the varsite, it is a done path, add to w
                    if (u[j].npath.back()->id == post_site_id)
                    {
                        w.push_back(u[j]);
                    } else if ((uint)std::accumulate(u[j].kmers_on_path.begin(), u[j].kmers_on_path.end(), 0) < 4*num_minis or
                               u[j].has_at_least_n_hit_minis_on_path(kmer_path_hit_counts[dir], num_minis)) {
                        // otherwise look up the last node in the max_path_index
                        // if it is there, then we have a path to extend by, 
                        // and extend the path for each of the maximal paths through the sub_varsite
                        // and add the updated path/prob pairs to v
                        map<uint32_t, vector<MaxPath>>::iterator it=bubble_paths.find(u[j].npath.back()->id);
                        if (it != bubble_paths.end())
                        {
                            //extend node path with the max paths seen before from this point
			    for (uint32_t i = 0; i != it->second.size(); ++i)
			    {
				v.push_back(u[j]);
			        v.back().extend(it->second[i]);
				update_kmers_on_node_path(v.back(), kmer_path_probs[dir]);
			    }
			} else {
			    //otherwise extend with the outnode of the last node in node_path
			    assert(u[j].npath.back()->outNodes.size() == 1); // if the back is the start of a bubble, should already be in index!
			    u[j].npath.push_back(u[j].npath.back()->outNodes[0]);
			    update_kmers_on_node_path(u[j], kmer_path_probs[dir]);
			    v.push_back(u[j]);
			}
                    } else {
			cout << "reject path " << u[j] << " which has " << std::accumulate(u[j].kmers_on_path.begin(), u[j].kmers_on_path.end(), 0);
                        cout << " minis on path and (at least " << num_minis << " have been hit) is ";
                        cout << u[j].has_at_least_n_hit_minis_on_path(kmer_path_hit_counts[dir], num_minis) << endl;
		    }
                }
                // once done for all of what was in u, set u = v
                u = v;
		v.clear();
            }

	    //assert(w.size()>0); // we have to have found some paths
	    assert(u.size()==0);

	    // at this point, all the paths should have up to date vectors of bools representing which kmers overlapped, 
	    // including the pre_site_id node and the post_site_id
	    // can work out which have enough minis which have hits

            for (uint n = 0; n!=w.size(); ++n)
            {
                if ((uint)std::accumulate(w[n].kmers_on_path.begin(), w[n].kmers_on_path.end(), 0) < 4*num_minis or 
		    w[n].has_at_least_n_hit_minis_on_path(kmer_path_hit_counts[dir], num_minis))
                {
		    cout << "keep path " << w[n] << " which has " << std::accumulate(w[n].kmers_on_path.begin(), w[n].kmers_on_path.end(), 0);
                    cout << " minis on path and (at least " << num_minis << " have been hit) is ";
                    cout << w[n].has_at_least_n_hit_minis_on_path(kmer_path_hit_counts[dir], num_minis) << endl;
		    u.push_back(w[n]);
                } else {
		    cout << "reject path " << w[n] << " which has " << std::accumulate(w[n].kmers_on_path.begin(), w[n].kmers_on_path.end(), 0);
                    cout << " minis on path and (at least " << num_minis << " have been hit) is ";
                    cout << w[n].has_at_least_n_hit_minis_on_path(kmer_path_hit_counts[dir], num_minis) << endl;
		}
	    }

            bubble_paths[pre_site_id] = u;
	    cout << "added " << u.size() << " paths to bubble_paths[" << pre_site_id << "]" << endl;
	    w.clear();
	    u.clear();
        }
    }

    // to finish return the vector of probs for 0 in bubble_paths
    vector<float> ret_values;
    for (uint n = 0; n != bubble_paths[0].size(); ++n)
    {
	ret_values.push_back(bubble_paths[0][n].get_prob(kmer_path_probs[dir]));
    }
    return ret_values;
}

/*void LocalPRG::write_max_path_to_vcf(const string& filepath)
{
    assert(max_path_index[0].size()==3);
    assert(max_path_index[0][0].size()>0);
    vector<LocalNode*> ref = prg.top_path();
    vector<LocalNode*> query = max_path_index[0][0][0].npath;
}*/

void LocalPRG::write_kmer_max_paths_to_fasta(const string& filepath, float e_rate)
{
    ofstream handle;
    handle.open (filepath);
    for (uint dir=0; dir!=2; dir++)
    {
	cout << now() << "find kmer max paths for dir " << dir << endl;
        vector<KmerNode*> kmp;
	kmp.reserve(30);
	float ppath = kmer_prg.find_max_path(dir, e_rate, kmp);
	vector<LocalNode*> lmp = localnode_path_from_kmernode_path(kmp);

        handle << ">" << name << "." << dir << "\t P(data|sequence)=" << ppath << endl;
        for (uint j = 0; j!= lmp.size(); ++j)
        {
            handle << lmp[j]->seq;
        }
        handle << endl;
    }
    handle.close();
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
    }
    handle.close();
    return;
}

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

