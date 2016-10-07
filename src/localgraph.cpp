#include <iostream>
#include <cstring>
#include <map>
#include <fstream>
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
    map<uint32_t, LocalNode*>::iterator it=nodes.find(id);
    if(it==nodes.end())
    {
        LocalNode *n;
        n = new LocalNode(seq, pos, id);
        nodes[id] = n;
        //cout << "Added node " << id << endl;
    } else {
        cerr << NODE_EXISTS_ERROR << endl;
	cerr << id << " " << pos << endl;
	exit(-1);
    }
    return;
}
void LocalGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    map<uint32_t, LocalNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, LocalNode*>::iterator to_it=nodes.find(to);
    if((from_it!=nodes.end()) && (to_it!=nodes.end()))
    {
	LocalNode *f = (nodes.find(from)->second);
    	LocalNode *t = (nodes.find(to)->second);
    	f->outNodes.push_back(t);
    	//cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    } else {
	cerr << NODE_MISSING_ERROR << endl;
	cerr << from << " " << to << endl;
	exit(-1);
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

/*uint32_t query_node_containing(uint32_t i)
{
    map<uint32_t, LocalNode*>::iterator it=nodes.begin();
    while (i<it->second->pos.end and it!=nodes.end())
    {
	++it;
    }
    if (it->second->pos.start<=i and it->second->pos.end>i)
    {
	return i;
    } else
	return 
}

vector<Path> LocalGraph::get_paths_from_node(uint32_t i, uint32_t len)
{
    vector<Path> v;
    deque<Interval> d = {Interval(i,15)};
    Path p = Path(
    for 
    return v;
}*/

set<Path> LocalGraph::extend_path(Path p)
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
	    cout << q << endl;
	    s.insert(q);
	}
    }
    return s;	
}
