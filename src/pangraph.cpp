#include <iostream>
#include <cstring>
#include <map>
#include <vector>
#include <fstream>
#include "utils.h"
#include "pannode.h"
#include "pangraph.h"
#include <cassert>
#include "minihit.h"

using namespace std;

PanGraph::~PanGraph()
{
  for (auto c: nodes)
  {
    delete c.second;
  }
}
/*void PanGraph::add_node (const uint32_t prg_id, const uint32_t read_id) // make cluster and optional parameter, by defining two functions
{
    map<uint32_t, PanNode*>::iterator it=nodes.find(prg_id);
    if(it==nodes.end())
    {
        PanNode *n;
        n = new PanNode(prg_id);
        nodes[prg_id] = n;
        //cout << "Added node " << *n << endl;
        n->add_read(read_id);
    } else {
        //cout << "Node " << prg_id << " was already in graph" << endl;
        it->second->add_read(read_id);
    }
    return;
}*/

void PanGraph::add_node (const uint32_t prg_id, const uint32_t read_id, const set<MinimizerHit*, pComp>& cluster)
{
    //cout << "prg: " << prg_id << " read: " << read_id << " cluster_size: " << cluster.size() << endl;
    for (set<MinimizerHit*, pComp>::iterator it = cluster.begin(); it != cluster.end(); ++it)
    {
            //cout << "read id and cluster read ids: " << read_id << " " << (*it)->read_id << endl;
            assert(read_id == (*it)->read_id); // the hits should correspond to the read we are saying...
	    assert(prg_id == (*it)->prg_id);
    }
    map<uint32_t, PanNode*>::iterator it=nodes.find(prg_id);
    if(it==nodes.end())
    {
        PanNode *n;
        n = new PanNode(prg_id);
        nodes[prg_id] = n;
	n->add_read(read_id);
        n->add_hits(cluster);
	//cout << "Added node " << *n << " and " << n->foundHits.size() << " hits" << endl;
    } else {
        //cout << "Node " << prg_id << " was already in graph" << endl;
        it->second->add_read(read_id);
	it->second->add_hits(cluster);
	//cout << "Added hits to node " << *(it->second) << " giving a total of " << it->second->foundHits.size() << " hits" << endl;
    }
    return;
}
void PanGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    map<uint32_t, PanNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, PanNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) and (to_it!=nodes.end()));
    PanNode *f = (nodes.find(from)->second);
    PanNode *t = (nodes.find(to)->second);
    f->outNodes.push_back(t);
    //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
}

void PanGraph::write_gfa (const string& filepath)
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
    // false if have different numbers of nodes
    //cout << "Graph == comparison:" << endl;
    //cout << "numbers of nodes: " << nodes.size() << " " << y.nodes.size() << endl;
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;}

    // false if have different nodes
    for ( const auto c: nodes)
    {
        //cout << "for node " << c.first;
        // if node id doesn't exist 
        map<uint32_t, PanNode*>::const_iterator it=y.nodes.find(c.first);
        if(it==y.nodes.end()) {//cout << "node id doesn't exist" << endl; 
            return false;}
        //cout << " found corresponding y node " << it->first;
        // or node entries are different
        if (!(*c.second == *(it->second))) {//cout << "node id " << c.first << " exists but has different values" << endl; 
            return false;}
        //cout << " such that " << *c.second << " == " << *(it->second) << endl;
    }
    // otherwise is true
    return true;
}

std::ostream& operator<< (std::ostream & out, PanGraph const& m) {
    for(map<uint32_t, PanNode*>::const_iterator it=m.nodes.begin(); it!=m.nodes.end(); ++it)
    {
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            out << it->second->id << "->" << it->second->outNodes[j]->id << endl;
        }
    }
    out << endl;
    return out ;
}
