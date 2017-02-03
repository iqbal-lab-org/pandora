#include <iostream>
#include <cstring>
#include <map>
#include <fstream>
#include <cassert>
#include <vector>
#include "localnode.h"
#include "localgraph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

LocalGraph::LocalGraph()
{
    // reserve space in index
    index.reserve(10);
}

LocalGraph::~LocalGraph()
{
  for (auto c: nodes)
  {
    delete c.second;
  }
}

void LocalGraph::add_node (const uint32_t& id, const string& seq, const Interval& pos)
{
    assert(seq.length() == pos.length); 
    map<uint32_t, LocalNode*>::iterator it=nodes.find(id);
    if(it==nodes.end())
    {
        LocalNode *n;
        n = new LocalNode(seq, pos, id);
        nodes[id] = n;
        //cout << "Added node " << id << endl;
    } else {
    	assert((it->second->seq == seq) && (it->second->pos == pos));
    }
    return;
}
void LocalGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    map<uint32_t, LocalNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, LocalNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) && (to_it!=nodes.end()));
    if((from_it!=nodes.end()) && (to_it!=nodes.end()))
    {
	LocalNode *f = (nodes.find(from)->second);
    	LocalNode *t = (nodes.find(to)->second);
	assert(f->pos.end<=t->pos.start || assert_msg(f->pos.end << ">" << t->pos.start << " so cannot add edge from node " << *f << " to node " << *t));
    	f->outNodes.push_back(t);
    	//cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    }
    return;
}

void LocalGraph::add_varsite (const uint8_t level, const uint32_t pre_site_id, const uint32_t post_site_id)
{
    assert(pre_site_id <= post_site_id);
    while (level >= index.size())
    {
        vector<pair<uint32_t, uint32_t>> levelv;
	levelv.reserve(400);
	//levelv = {};
	index.insert(index.end(), 1, levelv);
    }
    index[level].push_back(make_pair(pre_site_id, post_site_id));
    return;
}

void LocalGraph::write_gfa (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for(map<uint32_t, LocalNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        handle << "S\t" << it->second->id << "\t";
        if (it->second->seq == "")
	{
	    handle << "*";
	} else {
	    handle << it->second->seq;
	}
	handle << "\tRC:i:" << it->second->covg << endl;
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            handle << "L\t" << it->second->id << "\t+\t" << it->second->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

vector<Path> LocalGraph::walk(const uint32_t& node_id, const uint32_t& pos, const uint32_t& len)
{
    //cout << "walking graph from node " << node_id << " pos " << pos << " for length " << len << endl;
    // walks from position pos in node node for length len bases
    assert(nodes[node_id]->pos.start <= pos && nodes[node_id]->pos.end >= pos); // if this fails, pos given lies on a different node
    vector<Path> return_paths, walk_paths;
    return_paths.reserve(20);
    walk_paths.reserve(20);
    Path p,p2;
    deque<Interval> d;

    //cout << "pos+len: " << pos+len << " nodes[node_id]->pos.end: " << nodes[node_id]->pos.end << endl;
    if (pos+len <= nodes[node_id]->pos.end)
    {
        d = {Interval(pos, pos+len)};
        p.initialize(d);
        //cout << "return path: " << p << endl;
        return_paths.push_back(p);
        //cout << "return_paths size: " << return_paths.size() << endl; 
        return return_paths;
    }
    uint32_t len_added = min(nodes[node_id]->pos.end - pos, len);

    //cout << "len: " << len << " len_added: " << len_added << endl;
    if (len_added < len)
    {
	for (vector<LocalNode*>::iterator it = nodes[node_id]->outNodes.begin(); it!= nodes[node_id]->outNodes.end(); ++it)
	{
	    //cout << "Following node: " << (*it)->id << " to add " << len-len_added << " more bases" << endl;
	    walk_paths = walk((*it)->id,(*it)->pos.start, len-len_added);
	    //cout << "walk paths size: " << walk_paths.size() << endl;
	    for (vector<Path>::iterator it2 = walk_paths.begin(); it2!= walk_paths.end(); ++it2)
	    {
		// Note, would have just added start interval to each item in walk_paths, but can't seem to force result of it2 to be non-const
		//cout << (*it2) << endl;
		p2.initialize((*it2).path);
		p2.add_start_interval(Interval(pos, nodes[node_id]->pos.end));
		//cout << "path: " << p2 << " p2.length: " << p2.length << endl;
    		if (p2.length == len) {
		    return_paths.push_back(p2);
		}
	    }
	}
    }
    return return_paths;
}

bool LocalGraph::operator == (const LocalGraph& y) const
{
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;}

    // false if have different nodes
    for ( const auto c: nodes)
    {
        // if node id doesn't exist 
        map<uint32_t, LocalNode*>::const_iterator it=y.nodes.find(c.first);
        if(it==y.nodes.end()) {//cout << "node id doesn't exist" << endl; 
            return false;}
        // or node entries are different
        if (!(*c.second == *(it->second))) {//cout << "node id " << c.first << " exists but has different values" << endl; 
            return false;}
    }
    // otherwise is true
    return true;
}

