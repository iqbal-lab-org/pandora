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
        cout << "Added node " << id << endl;
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
    	cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
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
