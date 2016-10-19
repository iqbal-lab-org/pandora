#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pannode.h"
#include "pangraph.h"
//#include "minimizerhit.h"

using namespace std;

PanGraph::~PanGraph()
{
  for (auto c: nodes)
  {
    delete c.second;
  }
}
void PanGraph::add_node (const uint32_t prg_id, const uint32_t read_id)
{
    map<uint32_t, PanNode*>::iterator it=nodes.find(prg_id);
    if(it==nodes.end())
    {
        PanNode *n;
        n = new PanNode(prg_id);
        nodes[prg_id] = n;
        cout << "Added node " << *n << endl;
	n->add_read(read_id);
        //n->add_hits(cluster);
    } else {
        cout << "Node " << prg_id << " was already in graph" << endl;
        it->second->add_read(read_id);
	//it->second->add_hits(cluster);
    }
    return;
}
void PanGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    map<uint32_t, PanNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, PanNode*>::iterator to_it=nodes.find(to);
    if((from_it!=nodes.end()) && (to_it!=nodes.end()))
    {
	PanNode *f = (nodes.find(from)->second);
    	PanNode *t = (nodes.find(to)->second);
    	f->outNodes.push_back(t);
    	cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    } else {
	cout << "One of the nodes was missing - could not add edge" << endl;
    }

}
void PanGraph::write_gfa (string filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        handle << "S\t" << it->second->id << "\t*" << endl; //\tRC:i:" << it->second->foundHits.size() * k << endl;
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            handle << "L\t" << it->second->id << "\t+\t" << it->second->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

bool PanGraph::operator == (const PanGraph& y) const {
    if (nodes.size() != y.nodes.size()) {return false;}
    for(map<uint32_t, PanNode*>::const_iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	// check node exists in y
        map<uint32_t, PanNode*>::const_iterator f=y.nodes.find(it->first);
	if (f==y.nodes.end()) {return false;}
	// and that each edge from node to other node is also in y
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            pointer_values_equal<PanNode> eq = { it->second->outNodes[j] };
            if ( find_if(f->second->outNodes.begin(), f->second->outNodes.end(), eq) == f->second->outNodes.end() )
	    {return false;}
        }
    }
    return true;
}

/*std::ostream& operator<< (std::ostream & out, PanGraph const& m) {
    out << ;
    return out ;
}*/
