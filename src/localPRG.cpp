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
#include "errormessages.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

LocalPRG::LocalPRG (uint32_t i, string n, string p): next_id(0), buff(" "), next_site(5), id(i), name(n), seq(p), num_hits(2, 0)
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

    //kmer_paths.reserve(50000);
    //kmer_path_clashes.reserve(50000);
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
    assert(s.length()==p.length()||assert_msg("sequence length " << s.length() << " is not equal to path length " << p.length()));
    return s;
}

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
            if ((it->end > n->second->pos.start and it->start < n->second->pos.end) or 
		(it->start == n->second->pos.start and it->end == n->second->pos.end) or 
		(it->start == n->second->pos.start and it->length == 0 and it == --(p.path.end()) and n!=prg.nodes.begin()))
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
	// add end interval
	end_ids = build_graph(v.back(), end_ids, current_level);
    }
    if (start_id == 0)
    {
	assert(end_ids.size()==1);
    }
    return end_ids;
}

vector<Path> LocalPRG::shift(Path p)
{
    // returns all paths of the same length which have been shifted by one position along prg graph
    Path q;
    q = p.subpath(1,p.length()-1);
    vector<LocalNode*> n;
    vector<Path> return_paths;
    deque<Path> short_paths = {q};
    vector<Path> k_paths;
    Interval i;
    bool non_terminus;

    // first find extensions of the path
    while (short_paths.size() > 0)
    {
	p = short_paths.front();
	n = nodes_along_path(p);
	short_paths.pop_front();
	
	// if we can extend within the same localnode, do
	if (p.end < n.back()->pos.end)
        {
            p.path.back().end += 1;
	    p.end += 1;
	    p.path.back().length += 1;
	    k_paths.push_back(p);    
	} else if (p.end != (--(prg.nodes.end()))->second->pos.end) {
	    for (uint i=0; i!=n.back()->outNodes.size(); ++i)
	    {
		short_paths.push_back(p);
		short_paths.back().add_end_interval(Interval(n.back()->outNodes[i]->pos.start, n.back()->outNodes[i]->pos.start));
	    }
	}
    }

    // now check if by adding null nodes we reach the end of the prg
    for (uint i=0; i!= k_paths.size(); ++i)
    {
	short_paths = {k_paths[i]};
	non_terminus = false; // assume there all extensions are terminal i.e. reach end or prg

	while (short_paths.size() > 0)
	{
	    p = short_paths.front();
            n = nodes_along_path(p);
            short_paths.pop_front();
	    
	    if (n.back()->pos.end==(--(prg.nodes.end()))->second->pos.end)
            {
                return_paths.push_back(p);
	    } else if (n.back()->pos.end==p.end) {
                for (uint i=0; i!=n.back()->outNodes.size(); ++i)
		{
		    if (n.back()->outNodes[i]->pos.length == 0)
		    {
                        short_paths.push_back(p);
                        short_paths.back().add_end_interval(n.back()->outNodes[i]->pos);
		    } else {
			non_terminus = true;
		    }
		}
            } else {
		non_terminus = true;
	    }
	}
	if (non_terminus == true)
	{
	    return_paths.push_back(k_paths[i]);
	}
    }

    return return_paths;
}

void LocalPRG::minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k)
{
    cout << now() << "Sketch PRG " << name << " which has " << prg.nodes.size() << " nodes" << endl;

    // clean up after any previous runs
    // although note we can't clear the index because it is also added to by other LocalPRGs
    kmer_prg.clear();

    // declare variables
    vector<Path> walk_paths, shift_paths, v;
    deque<KmerNode*> current_leaves, end_leaves;
    deque<vector<Path>> shifts;
    deque<Interval> d;
    Path kmer_path;
    string kmer;
    uint64_t smallest;
    pair<uint64_t, uint64_t> kh;
    KmerHash hash;
    uint num_kmers_added = 0;
    KmerNode *kn, *new_kn;
    vector<KmerNode*>::iterator found;
    vector<LocalNode*> n;
    bool mini_found_in_window;

    // create a null start node in the kmer graph
    d = {Interval(0,0)};
    kmer_path.initialize(d);
    kn = kmer_prg.add_node(kmer_path);
    num_kmers_added += 1;

    // if this is a null prg, return the null kmergraph
    if (prg.nodes.size() == 1 and prg.nodes[0]->pos.length < k)
    {
	return;
    }

    /*// force every start kmer to be in the graph
    walk_paths = prg.walk(prg.nodes.begin()->second->id, 0, k);
    for (uint i=0; i!=walk_paths.size(); ++i)
    {
        kmer = string_along_path(walk_paths[i]);
        kh = hash.kmerhash(kmer, k);

	// add to index, kmer_prg and kmer_paths
	idx->add_record(min(kh.first, kh.second), id, walk_paths[i], (kh.first<kh.second));
	kn = kmer_prg.add_node_with_kh(walk_paths[i], min(kh.first, kh.second));
	current_leaves.push_back(kn);
	num_kmers_added += 1;
        kmer_prg.add_edge(kmer_prg.nodes[0]->path, walk_paths[i]);
    }*/

    // find first w,k minimizers
    walk_paths = prg.walk(prg.nodes.begin()->second->id, 0, w+k-1);
    for (uint i=0; i!=walk_paths.size(); ++i)
    {
        // find minimizer for this path 
        smallest = std::numeric_limits<uint64_t>::max();
	mini_found_in_window = false;
        for (uint32_t j = 0; j != w; j++)
        {
            kmer_path = walk_paths[i].subpath(j,k);
            if (kmer_path.path.size() > 0)
            {
                kmer = string_along_path(kmer_path);
                kh = hash.kmerhash(kmer, k);
                smallest = min(smallest, min(kh.first, kh.second));
            }
        }
        for (uint32_t j = 0; j != w; j++)
        {
	    kmer_path = walk_paths[i].subpath(j,k);
            if (kmer_path.path.size() > 0)
            {       
                kmer = string_along_path(kmer_path);
                kh = hash.kmerhash(kmer, k);
		n = nodes_along_path(kmer_path);
	
                if (prg.walk(n.back()->id, n.back()->pos.end, w+k-1).size() == 0)
                {
                    while (kmer_path.end >= n.back()->pos.end and n.back()->outNodes.size() == 1 and n.back()->outNodes[0]->pos.length == 0)
                    {
                        kmer_path.add_end_interval(n.back()->outNodes[0]->pos);
                        n.push_back(n.back()->outNodes[0]);
                    }
                }
	
                if (kh.first == smallest or kh.second == smallest)
                {
		    found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(kmer_path));
                    if (found == kmer_prg.nodes.end())
                    {
		        // add to index, kmer_prg
		        //cout << "add first minikmer for i:" << i << " j: " << j << " kmer: " << kmer << " kh:" << min(kh.first, kh.second)  << " and path: " << kmer_path << endl;
		        idx->add_record(min(kh.first, kh.second), id, kmer_path, (kh.first<=kh.second));
		        kn = kmer_prg.add_node_with_kh(kmer_path, min(kh.first, kh.second));	
		        num_kmers_added += 1;
		        if (mini_found_in_window == false)
		        {
                            kmer_prg.add_edge(kmer_prg.nodes[0]->path, kmer_path);
		        }
		        mini_found_in_window = true;
		        current_leaves.push_back(kn);
		    }
		}
	    }
	}
    }

    // while we have intermediate leaves of the kmergraph, for each in turn, explore the neighbourhood
    // in the prg to find the next minikmers as you walk the prg
    while (current_leaves.size()>0)
    {
	kn = current_leaves.front();
	current_leaves.pop_front();
	assert(kn->khash < std::numeric_limits<uint64_t>::max());
	//cout << "looking for outnodes of " << kn->path << endl;

        // find all paths which are this kmernode shifted by one place along the graph
	shift_paths = shift(kn->path);
	if (shift_paths.size() == 0)
	{
	    //assert(kn->path.start == 0); not true for a too short test, would be true if all paths long enough to have at least 2 minikmers on...
	    end_leaves.push_back(kn);
	}
	for (uint i=0; i!=shift_paths.size(); ++i)
	{
	    v = {shift_paths[i]};
	    shifts.push_back(v);
	}
	shift_paths.clear();

	while(shifts.size() > 0)
	{
	    v = shifts.front();
	    shifts.pop_front();	    
	    assert(v.back().length() == k);
	    kmer = string_along_path(v.back());
            kh = hash.kmerhash(kmer, k);
	    /*if (v.back().end == (--(prg.nodes.end()))->second->pos.end)
            {
		// this kmer is a terminus
		found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(v.back()));
                if (found == kmer_prg.nodes.end())
                {
		    idx->add_record(min(kh.first, kh.second), id, v.back(), (kh.first<=kh.second));
		    new_kn = kmer_prg.add_node(v.back());
	  	    kmer_prg.add_edge(kn, new_kn);
		    end_leaves.push_back(new_kn);
		    num_kmers_added += 1;
		} else {
		    kmer_prg.add_edge(kn, *found);
		}
	    } else */
	    //if (kh.first <= kn->khash or kh.second <= kn->khash) {
	    if (min(kh.first, kh.second) <= kn->khash) {
		// found next minimizer
		found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(v.back()));
		if (found == kmer_prg.nodes.end())
		{
		    //cout << "add minikmer: " << kmer << " kh:" << min(kh.first, kh.second)  << " and path: " << v.back() << " as smaller than " << kn->khash << endl;
		    idx->add_record(min(kh.first, kh.second), id, v.back(), (kh.first<=kh.second));
                    new_kn = kmer_prg.add_node_with_kh(v.back(), min(kh.first, kh.second));
                    kmer_prg.add_edge(kn, new_kn);
		    if (v.back().end == (--(prg.nodes.end()))->second->pos.end)
		    {
			end_leaves.push_back(new_kn);
		    } else {
                        current_leaves.push_back(new_kn);
		    }
                    num_kmers_added += 1;
		} else {
		    kmer_prg.add_edge(kn, *found);
		}
            } else if (v.size() == w) {
		// the old minimizer has dropped out the window, minimizer the w new kmers
		smallest = std::numeric_limits<uint64_t>::max();
		mini_found_in_window = false;
		//cout << "vsize is w" << endl;
                for (uint j = 0; j != w; j++)
                {
                    kmer = string_along_path(v[j]);
                    kh = hash.kmerhash(kmer, k);
                    smallest = min(smallest, min(kh.first, kh.second));
		    //cout << min(kh.first, kh.second) << " ";
                }
		for (uint j = 0; j != w; j++)
                {
		    kmer = string_along_path(v[j]);
                    kh = hash.kmerhash(kmer, k);
		    if (kh.first == smallest or kh.second == smallest)
		    {
			    found = find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(v[j]));
			    if (found == kmer_prg.nodes.end())
			    {
				//cout << "add minikmer for j: " << j << " kmer: " << kmer << " kh:" << min(kh.first, kh.second)  << " and path: " << v[j] << " since v.size()==" << v.size() << endl;
			        idx->add_record(min(kh.first, kh.second), id, v[j], (kh.first<=kh.second));
                                new_kn = kmer_prg.add_node_with_kh(v[j], min(kh.first, kh.second));

				// if there is more than one mini in the window, edge should go to the first, and from the first to the second
				if (mini_found_in_window == false)
                                {
                                    kmer_prg.add_edge(kn, new_kn);
                                }
                                mini_found_in_window = true;

				if (v.back().end == (--(prg.nodes.end()))->second->pos.end)
                                {
                                    end_leaves.push_back(new_kn);
                                } else {
                                    current_leaves.push_back(new_kn);
                                }
                                num_kmers_added += 1;
			    } else {
				if (mini_found_in_window == false)
				{
				    kmer_prg.add_edge(kn, *found);
				}
				mini_found_in_window = true;
			    }
		    }
		}
	    } else if (v.back().end == (--(prg.nodes.end()))->second->pos.end) {
		//cout << "reached terminus" << endl;
                end_leaves.push_back(kn);
	    } else {
		//cout << "not found a new mini or terminus so shift again";
		shift_paths = shift(v.back());
                for (uint i=0; i!=shift_paths.size(); ++i)
                {
		    shifts.push_back(v);
		    shifts.back().push_back(shift_paths[i]);
                }
                shift_paths.clear();
		//cout << " - shifts now has size " << shifts.size() << endl;
	    }
	}
    }

    // create a null end node, and for each end leaf add an edge to this terminus
    assert(end_leaves.size()>0);
    d = {Interval((--(prg.nodes.end()))->second->pos.end,(--(prg.nodes.end()))->second->pos.end)};
    kmer_path.initialize(d);
    kn = kmer_prg.add_node(kmer_path);
    num_kmers_added += 1;
    for (uint i=0; i!=end_leaves.size(); ++i)
    {
        kmer_prg.add_edge(end_leaves[i], kn);
    }

    // print, check and return
    //cout << "kmer prg: " << endl << kmer_prg << endl;
    //cout << "sort kmernodes topologically" << endl;
    kmer_prg.sort_topologically();
    //cout << "check prg" << endl;
    kmer_prg.check(num_kmers_added);
    return;
}

vector<KmerNode*> LocalPRG::find_kmernodes_on_localnode_path(vector<LocalNode*>& npath)
{

    deque<Interval>::const_iterator it;
    vector<LocalNode*>::const_iterator node;

    bool found, reject;
    vector<KmerNode*> nums;
    nums.reserve(500);

    // find the set of kmer_paths which overlap this node_path
    for (uint32_t n=1; n!=kmer_prg.nodes.size() - 1; ++n)
    {
        if (kmer_prg.nodes[n]->path.end >= npath[0]->pos.start and kmer_prg.nodes[n]->path.start <= npath.back()->pos.end) {
	    
            found = false;
            reject = false;

            it=kmer_prg.nodes[n]->path.path.begin();
            node=npath.begin();

            while (reject == false and found == false and node!=npath.end())
            {
                if (it!=kmer_prg.nodes[n]->path.path.end() and
                   ((it->end > (*node)->pos.start and it->start < (*node)->pos.end) or
                   (*it == (*node)->pos)))
                {
                    // then this node overlaps this interval of the kmer
                    // it needs to continue to overlap to end of kmer or node_path
                    found = true;
                    node++;
                    it++;
                    while (reject == false and node!=npath.end() and it!=kmer_prg.nodes[n]->path.path.end())
                    {
                        if (it!=kmer_prg.nodes[n]->path.path.end() and
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
                    if (node==npath.begin() and it!=kmer_prg.nodes[n]->path.path.end())
                    {
                        it++;
                    } else {
                        it=kmer_prg.nodes[n]->path.path.begin();
                        node++;
                    }
                }
            }
            // now if it was found and not rejected, add to nums;
            if (found == true and reject == false)
            {
		if (nums.size() > 0 and nums.back()->path.end == kmer_prg.nodes[n]->path.end)
		{
		    if (nums.back()->path.start > kmer_prg.nodes[n]->path.start)
		    {
			nums.pop_back();
		    }
		} else if (nums.size() > 0 and nums.back()->path.start == kmer_prg.nodes[n]->path.start)
                {
                    if (nums.back()->path.end < kmer_prg.nodes[n]->path.end)
                    {
                        nums.pop_back();
                    }
                } else {
		/*if (equal_except_null_nodes(nums.back()->path, kmer_prg.nodes[n]->path))
		{
		    if (nums.back()->path.start > kmer_prg.nodes[n]->path.start)
		    {
			nums.pop_back();
		    }
		} else {*/
                    nums.push_back(kmer_prg.nodes[n]);
		}
                //cout << "found kmer match for " << kmer_paths[n] << endl;
            }
        }
    }
    return nums;
}

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
    }

    cout << endl << "so localnode path found is: " << endl;
    for(uint i=0; i!=localnode_path.size(); ++i)
    {
	cout << *localnode_path[i] << " -> ";
    }
    cout << endl;
    return localnode_path;
}

void LocalPRG::update_covg_with_hit(MinimizerHit* mh)
{
    bool added = false;
    // update the covg in the kmer_prg
    for (uint i=0; i!=kmer_prg.nodes.size(); ++i)
    {
	if (kmer_prg.nodes[i]->path == mh->prg_path)
	{
	    kmer_prg.nodes[i]->covg[mh->strand] += 1;
	    cout << "kmernode " << i << " with path " << kmer_prg.nodes[i]->path << " now has " << mh->strand << " coverage " << kmer_prg.nodes[i]->covg[mh->strand] << endl;
	    added = true;
	    break;
	}
    }
    num_hits[mh->strand] += 1;
    assert(added == true || assert_msg("could not find kmernode corresponding to " << *mh));
}

void LocalPRG::write_kmer_max_paths_to_fasta(const string& filepath, float e_rate)
{
    ofstream handle;
    handle.open (filepath);

    vector<KmerNode*> kmp;
    kmp.reserve(800);
    float ppath;
    vector<LocalNode*> lmp;
    
    ppath = kmer_prg.find_max_path(e_rate, kmp);
    lmp = localnode_path_from_kmernode_path(kmp);
    handle << ">" << name << "\tlog P(data|sequence)=" << ppath  << endl;
    for (uint j = 0; j!= lmp.size(); ++j)
    {
        handle << lmp[j]->seq;
    }
    handle << endl;

    /*cout << now() << "find kmer max paths for reverse complement direction" << endl;
    kmp_b.clear();
    ppath_b = kmer_prg.find_max_path(0, e_rate, kmp_b);
    lmp = localnode_path_from_kmernode_path(kmp_b);
    handle << ">" << name << ".rc" << "\tlog P(data|sequence)=" << ppath_b << endl;
    for (uint j = lmp.size(); j!= 0; --j)
    {
        handle << rev_complement(lmp[j-1]->seq);
    }
    handle << endl;*/

    handle.close();
    return;
}

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

