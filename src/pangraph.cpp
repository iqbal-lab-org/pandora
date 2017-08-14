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

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

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
void PanGraph::add_edge (const uint32_t& from, const uint32_t& to, const uint& orientation, const uint& read_id)
{
    map<uint32_t, PanNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, PanNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) and (to_it!=nodes.end()));
    assert((orientation < (uint)4) || assert_msg("tried to add an edge with a rubbish orientation " << orientation << " which should be < 4"));
    PanNode *f = (nodes.find(from)->second);
    PanNode *t = (nodes.find(to)->second);
    //f->outNodes.push_back(t);

    assert(f->outNodes.size() == 4);
    assert(t->outNodes.size() == 4);
    assert(f->outNodeCounts.size() == 4);
    assert(t->outNodeCounts.size() == 4);

    //cout << orientation << " " << 3-orientation << endl;
    if (std::find(f->outNodes[orientation].begin(), f->outNodes[orientation].end(), t) == f->outNodes[orientation].end())
    {
	//cout << "Added edge (" << f->id << ", " << t->id << " " << orientation << ")" << endl;
	f->outNodes[orientation].push_back(t);
	f->outNodeCounts[orientation][t->id].push_back(read_id);
    } else {
	//cout << "Added count for (" << f->id << ", " << t->id << " " << orientation << ")" << endl;
	f->outNodeCounts[orientation][t->id].push_back(read_id);
    }

    if (std::find(t->outNodes[3-orientation].begin(), t->outNodes[3-orientation].end(), f) == t->outNodes[3-orientation].end())
    {
	//cout << "Added edge (" << t->id << ", " << f->id << " " << 3-orientation << ")" << endl;
        t->outNodes[3-orientation].push_back(f);
        t->outNodeCounts[3-orientation][f->id].push_back(read_id);
    } else {
	//cout << "Added count for (" << t->id << ", " << f->id << " " << 3-orientation << ")" << endl;
        t->outNodeCounts[3-orientation][f->id].push_back(read_id);
    }
    //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
}

void PanGraph::delete_edge(PanNode* from, PanNode* to, const uint& orientation)
{
    cout << "delete edge " << to->name << "->" << from->name << " " << 3-orientation << endl;
    to->outNodes[3-orientation].erase(std::remove(to->outNodes[3-orientation].begin(),
                                                  to->outNodes[3-orientation].end(),
                                                  from),
                                      to->outNodes[3-orientation].end());
    cout << "delete edge " << from->name << "->" << to->name << " " << orientation << endl;
    from->outNodes[orientation].erase(std::remove(from->outNodes[orientation].begin(),
                                                  from->outNodes[orientation].end(),
                                                  to),
                                      from->outNodes[orientation].end());

    cout << "replace with edges joining to next node in read" << endl;
    // for each compatible orientation of an edge B->C
    uint k = 0;
    if (orientation > 1)
    {
	k = 1;
    }
    while (k < 4)
    {
	// iterate over the map of outnotes of 'to'
	for (std::unordered_map<uint32_t,vector<uint32_t>>::iterator it=to->outNodeCounts[k].begin(); it!=to->outNodeCounts[k].end(); ++it)
	{
	    // for each read id on deleted edge, search for the read id for an edge from 'to' to outnode
	    for (vector<uint32_t>::iterator id=from->outNodeCounts[orientation][to->id].begin(); id!=from->outNodeCounts[orientation][to->id].end(); ++id)
	    {
		if (find(it->second.begin(), it->second.end(), *id) != it->second.end())
		{
		    // if we find one, add a new edge
		    cout << "added edge from " << from->id << "->" << it->first << ", " << orientation + k -3*(orientation>1) << " for read " << *id << endl;
		    add_edge(from->id, it->first, orientation + k -3*(orientation>1), *id);
		}
	    }
	}
	k += 2;
    }
    to->outNodeCounts[3-orientation].erase(from->id);
    from->outNodeCounts[orientation].erase(to->id);
}

void PanGraph::clean(const uint32_t& covg)
{
    // first check
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {  
        for (uint k=0; k<4; ++k)
        {
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
	    {
		if (std::find(it->second->outNodes[k][j-1]->outNodes[3-k].begin(), it->second->outNodes[k][j-1]->outNodes[3-k].end(), it->second) == it->second->outNodes[k][j-1]->outNodes[3-k].end())
		{
		    cout << "asymmetric edge between " << it->second->name << " and " << it->second->outNodes[k][j-1]->name << endl;
		}
	    }
	}
    }
        
    uint thresh = 0.025*covg;
    cout << now() << "Cleaning with threshold " << thresh << endl;
    bool found;

    // remove all edges with less than 2.5% coverage (mostly we want to remove all the singley occuring edges)
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	//cout << *(it->second) << endl;
	for (uint k=0; k<4; ++k)
	{
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
            {
	        //cout << "outnode " << it->second->outNodes[k][j-1]->id;
	        //cout << " has count " << it->second->outNodeCounts[it->second->outNodes[k][j-1]->id] << endl;
	        if (it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].size() < thresh)
	        {
		    // delete the edge
		    /*cout << "delete edge " << it->second->outNodes[k][j-1]->name << "->" << it->second->name << " " << 3-k << endl;
                    it->second->outNodes[k][j-1]->outNodes[3-k].erase(std::remove(it->second->outNodes[k][j-1]->outNodes[3-k].begin(),
                                                                                  it->second->outNodes[k][j-1]->outNodes[3-k].end(),
                                                                                  it->second),
                                                                      it->second->outNodes[k][j-1]->outNodes[3-k].end());
                    cout << "delete edge " << it->second->name << "->" << it->second->outNodes[k][j-1]->name << " " << k << endl;
		    it->second->outNodes[k].erase(it->second->outNodes[k].begin() + j-1);*/
		    delete_edge(it->second, it->second->outNodes[k][j-1], k);
	        }
            }
	}
    }
    // and again for higher coverage
    thresh = 0.05*covg;
    cout << now() << "Cleaning with threshold " << thresh << endl;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        for (uint k=0; k<4; ++k)
        {
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
            {
                if (it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].size() < thresh)
                {
                    delete_edge(it->second, it->second->outNodes[k][j-1], k);
                }
            }
        }
    }
    // remove unconnected nodes
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end();)
    {
	found = false;
	for (uint k=0; k<4; ++k)
	{
	    if (it->second->outNodes[k].size() > 0)
	    {
	        // found outnodes, or innodes
		found = true;
		break;
	    }
	}
	if (found == false)
	{
	    cout << "delete node " << it->second->name;
	    delete it->second;
	    nodes.erase(it++);
	    cout << " so pangraph now has " << nodes.size() << " nodes" << endl;
	} else {
	    ++it;
	}
    }	    

    // check again
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        for (uint k=0; k<4; ++k)
        {
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
            {
                if (std::find(it->second->outNodes[k][j-1]->outNodes[3-k].begin(), it->second->outNodes[k][j-1]->outNodes[3-k].end(), it->second) == it->second->outNodes[k][j-1]->outNodes[3-k].end())
                {
                    cout << "asymmetric edge between " << it->second->name << " and " << it->second->outNodes[k][j-1]->name << endl;
                }
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
	cout << it->second->id << endl;
	
        handle << "S\t" << it->second->name << "\t*" << endl; //\tRC:i:" << it->second->foundHits.size() * k << endl;
        for (uint32_t j=0; j<it->second->outNodes[3].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[3][j]->name << "\t+\t0M\tRC:i:" << it->second->outNodeCounts[3][it->second->outNodes[3][j]->id].size() << endl;
        }
	for (uint32_t j=0; j<it->second->outNodes[1].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[1][j]->name << "\t-\t0M\tRC:i:" << it->second->outNodeCounts[1][it->second->outNodes[1][j]->id].size() << endl;
        }
    }
    handle.close();
}

std::ostream& operator<< (std::ostream & out, PanGraph const& m) {
    for(map<uint32_t, PanNode*>::const_iterator it=m.nodes.begin(); it!=m.nodes.end(); ++it)
    {
        for (uint32_t j=0; j<it->second->outNodes[3].size(); ++j)
        {
            out << "+" << it->second->id << " -> +" << it->second->outNodes[3][j]->id << endl;
        }
	for (uint32_t j=0; j<it->second->outNodes[1].size(); ++j)
        {
            out << "+" << it->second->id << " -> -" << it->second->outNodes[1][j]->id << endl;
        }
	
    }
    out << endl;
    return out ;
}
