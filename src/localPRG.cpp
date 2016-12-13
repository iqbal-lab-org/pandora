#include <iostream>
#include <fstream>
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
#include "path.h"
#include "minihits.h"
#include "minihit.h"
#include "inthash.h"
#include "pannode.h"
#include "utils.h"
#include "errormessages.h"

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
        // find the appropriate node of the prg
        for (map<uint32_t, LocalNode*>::const_iterator n=prg.nodes.begin(); n!=prg.nodes.end(); ++n)
        {
            if ((it->end > n->second->pos.start and it->start < n->second->pos.end) or (it->start == n->second->pos.start and it->end == n->second->pos.end))
            {
		v.push_back(n->second);
            } else if (it->end < n->second->pos.start)
            { break;} // because the local nodes are labelled in order of occurance in the linear prg string, we don't need to search after this
        }
    }
    return v;
}

vector<Interval> LocalPRG::splitBySite(const Interval& i)
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
    //assert(v.size()==1 or v.size()==3);

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
	    vector<uint32_t> w = build_graph(v[j], mid_ids, current_level+1);
	    end_ids.insert(end_ids.end(), w.begin(), w.end());
	}
	// add this varsite to index
	prg.add_varsite (current_level + 1, pre_site_id, next_id);
	// add end interval
	end_ids = build_graph(v.back(), end_ids);
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
            //cout << "Node has been deemed long enough and " << it->second->pos.start << " <= " << i << " <= " << it->second->pos.end - w - k + 1 << endl; 
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
	    if (i > it->second->pos.end+1-w-k)
	    {break;} // too close to the end of this node for this method
	}

        for (uint32_t i=it->second->sketch_next; i!=it->second->pos.end;)
        {
            walk_paths = prg.walk(it->second->id, i, w+k-1);
            cout << "for id, i: " << it->second->id << ", " << i << " found " << walk_paths.size() << " paths" << endl;
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
            if (i > it->second->pos.end+1-w-k)
            {break;} // too close to the end of this node for this method
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

void LocalPRG::infer_most_likely_prg_paths_for_corresponding_pannode(const PanNode* pnode, uint32_t k, float e_rate)
{
    // note that within the foundHits for a pnode, hits which map to same read, prg and strand will be grouped together
    // we want to know for each prg minimizer, how much support there is from the reads
    // start by counting how many hits against each minimizing_kmer of prg
    cout << "start by counting how many hits against each minimizing kmer in prg" << endl;
    vector<uint32_t> kmer_path_hit_counts(kmer_paths.size(),0);
    cout << "there are " << kmer_paths.size() << " minimizing kmers to count hits against" << endl;
    set<MinimizerHit*, pComp_path>::iterator mh_previous = pnode->foundHits.begin();
    uint32_t num_hits = 1, num_kmers=0;
    for (set<MinimizerHit*, pComp_path>::iterator mh = (pnode->foundHits.begin()); mh != pnode->foundHits.end();)
    {
        mh++;
        if (mh==pnode->foundHits.end() or (*mh)->strand != (*mh_previous)->strand or !((*mh)->prg_path == (*mh_previous)->prg_path))
        {
            for (uint32_t i=num_kmers; i!=kmer_paths.size(); ++i)
            {
                if (kmer_paths[i]==(*mh_previous)->prg_path)
                {
                    //cout << "found path " << kmer_paths[i] << " == " << (*mh_previous)->prg_path << endl;
                    kmer_path_hit_counts[i] = num_hits;
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
    uint32_t sum_of_elems = 0;
    for (uint32_t n : kmer_path_hit_counts)
    {sum_of_elems += n;}
    cout << "after adding counts of the " << pnode->foundHits.size() << " hits, have still got " << sum_of_elems << " hits added to tallies" << endl;
    assert(sum_of_elems==pnode->foundHits.size());

    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    cout << "next work out prob of seeing this number of hits against a kmer assuming it is truly present" << endl;
    float p_thresh=0.5;
    vector<float> kmer_path_probs;
    float p_kmer, p=1/exp(e_rate*k), p_max=0, p_min=1;
    uint32_t n = pnode->foundReads.size(), big_p_count=0;
    cout << "n: " << n << ", p: " << p << endl;
    cout << "count:0 " << nchoosek(n,0) << " * " << pow(p,0) << " * " << pow(1-p,n-0) << endl;
    cout << "count:1 " << nchoosek(n,1) << " * " << pow(p,1) << " * " << pow(1-p,n-1) << endl;
    for (uint32_t i=0; i!=kmer_path_hit_counts.size(); ++i)
    {
        p_kmer = nchoosek(n, kmer_path_hit_counts[i])*pow(p,kmer_path_hit_counts[i])*pow(1-p,n-kmer_path_hit_counts[i]);
        kmer_path_probs.push_back(p_kmer);
        p_max = max(p_max, p_kmer);
        p_min = min(p_min, p_kmer);
        if(p_kmer>p_thresh)
        {
            //cout << p_kmer << " ";
            big_p_count+=1;
        }
    }
    cout << endl;
    cout << "found " << big_p_count << " probs > " << p_thresh << " with max and min probs " << p_max << ", " << p_min << endl;

    //now we iterate through the graph from the outmost level to the lowest level working out the most likely path(s)
    //need 2 data structures, one to remember what the most probable path(s) were for var sites already considered
    //the first is max_path_index, stored by the LocalPRG class
    //and the second remembering which minimizing kmers have been seen/used before so not included multiple times in the probability
    vector<bool> kmer_considered_before(kmer_paths.size(),false);

    // start with the outmost level
    uint8_t max_level = 0;
    for (auto const& element : prg.index) {
        max_level = max(max_level, element.first);
    }
    // and for each level..
    for (uint level = max_level; level <= max_level; --level)
    {
	cout << "Looking at level " << level << endl;
        // ...for each varsite at this level...
        for (uint i = 0; i!=prg.index[level].size(); ++i)
        {
            // ...find the maximally probable paths through varsite
            vector<pair<vector<LocalNode*>, float>> u, v, w, z;
            vector<LocalNode*> x;
	    float p_new = 1;

            // add the first node of each alternate allele for the varsite to a vector
            uint32_t pre_site_id = prg.index[level][i].first, post_site_id = prg.index[level][i].second;
	    if (level == 0)
	    {
		assert(pre_site_id == 0);
		cout << "finally find best path through whole prg" << endl;
		x.push_back(prg.nodes[pre_site_id]);
		u.push_back(make_pair(x, 1));
		x.clear();
	    } else {
	        cout << "Looking at varsite number " << i << " at this level, from the outnodes of " << pre_site_id << " to the innodes of " << post_site_id << endl;
		// we want the index to be inclusive of first node/pre_site_id prob, but exclusive of end node prob
		x.push_back(prg.nodes[pre_site_id]);
		u.push_back(make_pair(x, 1));
		x.clear();
                /*for (uint j = 0; j!=prg.nodes[pre_site_id]->outNodes.size(); ++j)
                {
                    x.push_back(prg.nodes[pre_site_id]->outNodes[j]);
                    u.push_back(make_pair(x, 1));
                    x.clear();
                }*/
	    }

	    assert (u.size()>0); // we just added either start node, or branching outnodes of a varsite to it!

            // then until we reach the end varsite:
            while (u.size()>0)
            {
		cout << "restart loop with u.size() == " << u.size() << endl;
                // for each pair in u
                for (uint j = 0; j!=u.size(); ++j)
                {
                    // if the path in the pair ends at the end of the varsite, it is a done path, add to w
                    if (u[j].first.back()->id == post_site_id)
                    {
                        w.push_back(u[j]);
                    } else {
                        // otherwise look up the last node in the max_path_index
                        // if it is there, then we can multiply the running total prob for this allele, 
                        // and extend the path for each of the maximal paths through the sub_varsite
                        // and add the updated path/prob pairs to v
                        map<uint32_t, vector<pair<vector<LocalNode*>, float>>>::iterator it=max_path_index.find(u[j].first.back()->id);
                        if (it != max_path_index.end())
                        {
                            //extend node path with the max paths seen before from this point
                            //NOTE need to find prob of start node here toooooooo
                            cout << "found " << u[j].first.back()->id << " in max_path_index with first prob " << it->second[0].second << " and " << it->second.size() << " paths" << endl;
                            for (uint n = 0; n!=it->second.size(); ++n)
                            {
                                v.push_back(u[j]);
                                v.back().first.insert(v.back().first.end(), it->second[n].first.begin(), it->second[n].first.end());
                                v.back().second = v.back().second * it->second[n].second;
                            }
                        } else {
                            // if not, then work out which minhits overlap this node, and multiply their probabilities, 
                            // then add the outnodes to get a new path/prob pair to add to v
                            // during this process, update kmer_considered_before to reflect the minihits now used in probabilities
                            // note that although this node is 'at the end of a site', it may also be the start of another site and have multiple outNodes
                            // also save this all to the index for repeat calls
                            cout << "did not find " << u[j].first.back()->id << " in max_path_index, so work out prob of node from scratch" << endl;
			    
                            for (uint32_t n=0; n!=kmer_paths.size(); ++n)
			    {
				// if we've already used this kmer
				if(kmer_considered_before[n]==false)
				{
				    for (deque<Interval>::const_iterator it2=kmer_paths[n].path.begin(); it2!=kmer_paths[n].path.end(); ++it2)
				    {
				        if (it2->end > u[j].first.back()->pos.start and it2->start < u[j].first.back()->pos.end)
				        {
					    // then this node overlaps this kmer
					    u[j].second = u[j].second * kmer_path_probs[n];
					    p_new = p_new * kmer_path_probs[n];
					    kmer_considered_before[n] = true;
					    break;
				        }
				    }
				}
			    }
			    
			    x.push_back(u[j].first.back()); // for index

			    for (uint32_t n=0; n!=u[j].first.back()->outNodes.size(); ++n)
			    {
				v.push_back(u[j]);
                                v.back().first.push_back(u[j].first.back()->outNodes[n]);
			    }
			    for (uint32_t n=0; n!=u[j].first.back()->outNodes.size(); ++n)
			    {	
				// and also add to lookup index to avoid repeat calculations
				// but again, if prob very small, only add the first to avoid explosion of stored paths
				z.push_back(make_pair(x, p_new));
				z.back().first.push_back(u[j].first.back()->outNodes[n]);
				if (p_new > 0 && p_new < 0.01)
				{
				    break;
				}
			    }
			    max_path_index[u[j].first.back()->id] = z;
                            x.clear();
			    z.clear();
			    p_new = 1;
			    

                        }
                    }
                }
                // once done for all of what was in u, set u = v
                u = v;
		v.clear();
            }

	    assert(w.size()>0); // we have to have found some paths

	    // if we are at level 0, we will want to include the probability of the last node as well
	    float p_last = 1;
	    if (level == 0)
	    {
		for (uint32_t n=0; n!=kmer_paths.size(); ++n)
                {
                    // if we've not already used this kmer
                    if(kmer_considered_before[n]==false)
                    {
			// if it overlaps the last node (wlog last of the first path in w)
                    	for (deque<Interval>::const_iterator it2=kmer_paths[n].path.begin(); it2!=kmer_paths[n].path.end(); ++it2)
                        {
                            if (it2->end > w[0].first.back()->pos.start and it2->start < w[0].first.back()->pos.end)
                            {
                                // then this node overlaps this kmer
                                p_last = p_last * kmer_path_probs[n];
                                kmer_considered_before[n] = true;
                                break;
                            }
                        }
                    }
                }
	    }

            // when u empty, should have final set in w and can work out the max of the probs
            float max_prob = 0;
            // find max for all probs we could work out values for
            for (uint n = 0; n!=w.size(); ++n)
            {
                if (w[n].second != 1) // if no mini kmers in prg overlapping node, can't define prob data came from that node, 
                                                   // so set the prob to 1 intially, then set to max path value 
                                                   // only happens for really short paths
                {
                    max_prob = max(max_prob, w[n].second);
                }
            }
	    cout << "max_prob for paths at this varsite: " << max_prob << endl;

            // now add paths achieving max, or if none add paths cannot judge to max_path_index
	    // note that in the case pre_site_id == 0 and level == 0, we may overwrite a previous entry to the index
	    // if the probability has dropped really low, only add the first path, to avoid exploding memory stores
            max_path_index[pre_site_id] = u; // know u is empty vector
            for (uint n = 0; n!=w.size(); ++n)
            {
                if ((max_prob == 0 and w[n].second == 1) or w[n].second == max_prob)
                {
		    w[n].second = w[n].second * p_last; // note that p_last is 1 unless we are on the last one
                    max_path_index[pre_site_id].push_back(w[n]);
		    if (max_prob < 0.01)
		    {
			break;
		    }
                }
		
            }
	    assert(max_path_index[pre_site_id].size()>0);
        }
    }
    return;
}

void LocalPRG::write_max_paths_to_fasta(const string& filepath)
{
    assert(max_path_index.size()>0);
    map<uint32_t, vector<pair<vector<LocalNode*>, float>>>::iterator it=max_path_index.find(0);
    assert(it!=max_path_index.end());
    assert(max_path_index[0].size()>0);

    ofstream handle;
    handle.open (filepath);
    for (uint i = 0; i!= max_path_index[0].size(); ++i)
    {
        handle << ">" << id << "." << i << "\t P(data|sequence)=" << max_path_index[0][i].second << endl;
	for (uint j = 0; j!= max_path_index[0][i].first.size(); ++j)
	{
            handle << max_path_index[0][i].first[j]->seq;
        }
        handle << endl;
    }
    handle.close();
    return;
}
