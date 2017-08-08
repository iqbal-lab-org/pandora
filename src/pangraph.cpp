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

void PanGraph::add_node (const uint32_t prg_id, const string prg_name, const uint32_t read_id, const set<MinimizerHit*, pComp>& cluster)
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
        n = new PanNode(prg_id, prg_name);
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
    //f->outNodes.push_back(t);

    if (std::find(f->outNodes.begin(), f->outNodes.end(), t) == f->outNodes.end())
    {
	f->outNodes.push_back(t);
	f->outNodeCounts[t->id] = 1;
    } else {
	f->outNodeCounts[t->id] += 1;
    }
    //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
}

void PanGraph::clean(const uint32_t& covg)
{
    uint thresh = 0.025*covg;
    cout << now() << "Cleaning with threshold " << thresh << endl;
    bool found;

    // remove all edges with less than 2.5% coverage (mostly we want to remove all the singley occuring edges)
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	//cout << *(it->second) << endl;
        for (uint32_t j=it->second->outNodes.size(); j!=0; --j)
        {
	    //cout << "outnode " << it->second->outNodes[j-1]->id;
	    //cout << " has count " << it->second->outNodeCounts[it->second->outNodes[j-1]->id] << endl;
	    if (it->second->outNodeCounts[it->second->outNodes[j-1]->id] < thresh)
	    {
		// delete the edge
		//cout << "delete edge " << it->second->id << "->" << it->second->outNodes[j-1]->id << endl;
		it->second->outNodes.erase(it->second->outNodes.begin() + j-1);
	    }
        }
    }
    //cout << "now remove isolated nodes" << endl;
    // remove unconnected nodes
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	if (it->second->outNodes.size() == 0)
	{
	    // no outnodes, look for innodes
	    found = false;
	    for(map<uint32_t, PanNode*>::iterator it2=nodes.begin(); it2!=nodes.end(); ++it2)
	    {
		if (std::find(it2->second->outNodes.begin(), it2->second->outNodes.end(), it->second) != it2->second->outNodes.end())
    		{
		    found = true;
		    break;
	    	}
	    }
	    if (found == false)
	    {
	        //cout << "delete node " << it->second->id << endl;
	        delete it->second;
	        nodes.erase(it);
	    }
	}
    }	    
}

bool PanGraph::operator == (const PanGraph& y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {
        return false;
    }

    // false if have different nodes
    for ( const auto c: nodes)
    {
        // if node id doesn't exist 
        map<uint32_t, PanNode*>::const_iterator it=y.nodes.find(c.first);
        if(it==y.nodes.end()) {
            return false;
	}
        // or node entries are different
        if (!(*c.second == *(it->second))) {
            return false;
	}
    }
    // otherwise is true
    return true;
}

void PanGraph::write_gfa (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        handle << "S\t" << it->second->name << "\t*" << endl; //\tRC:i:" << it->second->foundHits.size() * k << endl;
        for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[j]->name << "\t+\t0M\tRC:i:" << it->second->outNodeCounts[it->second->outNodes[j]->id] << endl;
        }
    }
    handle.close();
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
