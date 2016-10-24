#include <iostream>
#include <cstring>
#include <map>
#include <fstream>
#include <cassert>
#include <vector>
#include "localnode.h"
#include "localgraph.h"
#include "errormessages.h"

using namespace std;

LocalGraph::~LocalGraph()
{
  for (auto c: nodes)
  {
    delete c.second;
  }
}

void LocalGraph::add_node (const uint32_t& id, const string& seq, Interval pos)
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
    	assert((it->second->seq == seq) && (it->second->pos == pos)); // NODE_EXISTS_ERROR
    }
    return;
}
void LocalGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    map<uint32_t, LocalNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, LocalNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) && (to_it!=nodes.end())); // NODE_MISSING_ERROR
    if((from_it!=nodes.end()) && (to_it!=nodes.end()))
    {
	LocalNode *f = (nodes.find(from)->second);
    	LocalNode *t = (nodes.find(to)->second);
    	f->outNodes.push_back(t);
    	//cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    }
}

void LocalGraph::write_gfa (string filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
    for(map<uint32_t, LocalNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        handle << "S\t" << it->second->id << "\t" << it->second->seq << "\tRC:i:" << it->second->covg << endl;
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            handle << "L\t" << it->second->id << "\t+\t" << it->second->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

vector<Path> LocalGraph::walk(uint32_t node_id, uint32_t pos, uint32_t len)
{
    //cout << "walking graph from node " << node_id << " pos " << pos << " for length " << len << endl;
    // walks from position pos in node node for length len bases
    vector<Path> return_paths, walk_paths;
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
		} else {
		    //cout << "PATH TOO SHORT " << p2 << endl;
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

/*uint32_t query_node_containing(uint32_t i)
{
    map<uint32_t, LocalNode*>::iterator it=nodes.begin();
    while (i<it->second->pos.end and it!=nodes.end())
    {
	++it;
    }
    if (it->second->pos.start<=i and it->second->pos.end>i)
    {
	return it->first;
    } else {
	cerr << NOT_IN_NODE_RANGE_ERROR << endl;
        exit(-1);
    }
}

vector<Path> LocalGraph::get_paths_from_position(uint32_t i, uint32_t len)
{
    vector<Path> v;
    deque<Interval> d = {Interval(i,15)};
    Path p = Path(
    for 
    return v;
}*/

/*set<Path> LocalGraph::extend_path(Path p)
{
    set<Path> s;
    // adds next interval(s) to end of path
    map<uint32_t, LocalNode*>::iterator it = nodes.begin();
    while (!((it->second)->pos==p.path.back()) and it!=nodes.end())
    {
        ++it;
    }

    if (it==nodes.end())
    {
	cerr << INTERVAL_MISSING_FROM_PATH_ERROR << endl;
        cerr << p.path.back() << endl;
        exit(-1);
    } else {
	for (uint32_t i=0; i!=it->second->outNodes.size(); ++i)
	{
	    Path q = p;
	    q.add_end_interval(it->second->outNodes[i]->pos);
	    //cout << q << endl;
	    s.insert(q);
	}
    }
    return s;	
}*/
