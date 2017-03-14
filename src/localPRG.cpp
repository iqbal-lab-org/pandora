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
    assert(s.length()==p.length);
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

void LocalPRG::minimizer_sketch (Index* idx, const uint32_t w, const uint32_t k)
{
    cout << now() << "Sketch PRG " << name << " which has " << prg.nodes.size() << " nodes" << endl;

    // clean up after any previous runs
    kmer_prg.clear();
    //kmer_paths.clear();
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
    uint num_kmers_added = 0;

    // create a null start node in the kmer graph
    d = {Interval(0,0)};
    kmer_path.initialize(d);
    prg.nodes.begin()->second->prev_kmer_paths.insert(kmer_path);
    kmer_prg.add_node(kmer_path);
    num_kmers_added += 1;

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
        num_kmers_added += 1;
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
		cout << "compare last kmer " << min(kh.first, kh.second) << " with smallest kmer " << smallest << endl;
                if(kh.first <= smallest or kh.second <= smallest)
                {
		    idx->add_record(min(kh.first, kh.second), id, kmer_path, (kh.first>kh.second));
		    num_kmers_added += 1;
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
	        cout << "found " << walk_paths.size() << " walk paths" << endl;
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
	
                            if ((kh.first == smallest or kh.second == smallest) and // or n.back() == (--(prg.nodes.end()))->second) and 
				(find_if(kmer_prg.nodes.begin(), kmer_prg.nodes.end(), condition(kmer_path))==kmer_prg.nodes.end()))
                            {
			        // add to index, kmer_prg and kmer_paths
			        idx->add_record(min(kh.first, kh.second), id, kmer_path, (kh.first>kh.second));
                                num_kmers_added += 1;
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
	    //cout << "add previously unfound end kmer" << endl;
            kmer = string_along_path(walk_paths[i]);
            kh = hash.kmerhash(kmer, k);

            // add to index, kmer_prg and kmer_paths
            idx->add_record(min(kh.first, kh.second), id, walk_paths[i], (kh.first>kh.second));
	    n = nodes_along_path(walk_paths[i]);
	    kmer_prg.add_node(walk_paths[i]);
	    set<Path> s = n[0]->prev_kmer_paths;
	    //cout << "prev_kmer_paths from n[0] " << *n[0] << " has size " << s.size() << endl;
	    set<Path> t, u;
	    //cout << "while back " << *(kmer_prg.nodes.back()) << " has no innodes" << endl;
	    while (kmer_prg.nodes.back()->inNodes.size() == 0 and s.size() > 0)
	    {
		//cout << "try to find an innode from a set of size " << s.size() << endl;
		for (set<Path>::iterator l = s.begin(); l!=s.end(); ++l)
		{
		    //cout << "consider adding edge to path " << *l << endl;
		    if ((*l < walk_paths[i]) and !(walk_paths[i].is_branching(*l)))
                    {
                        kmer_prg.add_edge(*l, walk_paths[i]);
		    } else if ((walk_paths[i] < *l) and (!(walk_paths[i].is_branching(*l)) or (kmer_prg.nodes.back()->inNodes.size()==0 and t.size()==0))) {
			//cout << "consider innodes instead" << endl;
			u = kmer_prg.get_innodes(*l);
			assert(u.size() > 0);
		        t.insert(u.begin(), u.end());
		    } else {
			//cout << walk_paths[i] << " and " << *l << " are not compatible" << endl;
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
            num_kmers_added += 1;
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
    num_kmers_added += 1;

    for (set<Path>::iterator l = ends.begin(); l!=ends.end(); ++l)
    {
        kmer_prg.add_edge(*l, kmer_path);
    }

    // print, check and return
    //cout << "kmer prg: " << endl << kmer_prg << endl;
    kmer_prg.check(num_kmers_added);
    return;
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
    return localnode_path;
}

void LocalPRG::update_covg_with_hit(MinimizerHit* mh)
{
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

std::ostream& operator<< (std::ostream & out, LocalPRG const& data) {
    out << data.name;
    return out ;
}

