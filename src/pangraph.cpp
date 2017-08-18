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

    if (std::find(f->outNodes[orientation].begin(), f->outNodes[orientation].end(), t) == f->outNodes[orientation].end())
    {
	//cout << "Added edge (" << f->id << ", " << t->id << " " << orientation << ")" << endl;
	f->outNodes[orientation].push_back(t);
	f->outNodeCounts[orientation][t->id].push_back(read_id);
    } else {
	//cout << "Added count for (" << f->id << ", " << t->id << " " << orientation << ")" << endl;
	f->outNodeCounts[orientation][t->id].push_back(read_id);
    }

    if (std::find(t->outNodes[rev_orient(orientation)].begin(), t->outNodes[rev_orient(orientation)].end(), f) == t->outNodes[rev_orient(orientation)].end())
    {
	//cout << "Added edge (" << t->id << ", " << f->id << " " << 3-orientation << ")" << endl;
        t->outNodes[rev_orient(orientation)].push_back(f);
        t->outNodeCounts[rev_orient(orientation)][f->id].push_back(read_id);
    } else {
	//cout << "Added count for (" << t->id << ", " << f->id << " " << 3-orientation << ")" << endl;
        t->outNodeCounts[rev_orient(orientation)][f->id].push_back(read_id);
    }
    //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
}

void PanGraph::delete_edge(PanNode* from, PanNode* to, const uint& orientation)
{
    cout << now() << "delete edge " << to->name << "->" << from->name << " " << rev_orient(orientation) << endl;
    to->outNodes[rev_orient(orientation)].erase(std::remove(to->outNodes[rev_orient(orientation)].begin(),
                                                  to->outNodes[rev_orient(orientation)].end(),
                                                  from),
                                      to->outNodes[rev_orient(orientation)].end());
    cout << now() << "delete edge " << from->name << "->" << to->name << " " << orientation << endl;
    from->outNodes[orientation].erase(std::remove(from->outNodes[orientation].begin(),
                                                  from->outNodes[orientation].end(),
                                                  to),
                                      from->outNodes[orientation].end());
    to->outNodeCounts[rev_orient(orientation)].erase(from->id);
    from->outNodeCounts[orientation].erase(to->id);
}

uint PanGraph::rev_orient(const uint& orientation)
{
    // 3 A  -> B  = B- -> A- 0
    // 2 A- -> B  = B- -> A  2
    // 0 A- -> B- = B  -> A  3
    // 1 A  -> B- = B  -> A- 1

    uint r_orientation = orientation;
    if (orientation == 0)
    {
        r_orientation = 3;
    } else if (orientation == 3)
    {
        r_orientation = 0;
    }
    return r_orientation;
}

void PanGraph::check_graph_symmetry()
{
    // checks that all out edges A->B have a corresponding in edge B<-A
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        for (uint k=0; k<4; ++k)
        {
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
            {
                assert(std::find(it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].begin(), it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].end(), it->second) != it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].end() || assert_msg("asymmetric edge between " << it->second->name << " and " << it->second->outNodes[k][j-1]->name));
            }
        }
    }
}

void PanGraph::remove_low_covg_nodes(const uint& thresh)
{
    // if A <-x- B -y-> C then the resulting equivalent A -> C has orientation
    //       y
    //    0 1 2 3
    //  0 . 1 . 3
    //  1 1 . 3 .
    //x 2 . 0 . 2 
    //  3 0 . 2 .
    //

    cout << now() << "Node cleaning with threshold " << thresh << endl;
    PanNode *from;
    PanNode *to;
    bool found_from, found_to;
    uint new_k;
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	if (it->second->foundReads.size() < thresh)
	{
	    // for each read, search for an in edge and an out edge corresponding to the read
	    for (vector<uint>::iterator id = it->second->foundReads.begin(); id!= it->second->foundReads.end(); ++id)
	    {
		found_from = false;
		found_to = false;
		new_k = 0;
	        for (uint k=0; k<4; ++k)
	        {
		    for (uint j=it->second->outNodes[k].size(); j!=0; --j)
		    {
			cout << "looking for read " << *id << endl;
		        cout << "orientation " << k << " outnode " << it->second->outNodes[k][j-1]->name << " is covered by reads ";
			for (uint i=0; i!=it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].size(); ++i)
			{
			    cout << it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id][i] << " ";
			}
			cout << endl;
			if (find(it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].begin(), it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].end(), *id) != it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].end())
			{
			    if (k == 0 or k == 2)
			    {
				cout << "found from edge " << it->second->name << " -> " << it->second->outNodes[k][j-1]->name << " " << k << " for read " << *id << endl;
			        found_from = true;
				from = it->second->outNodes[k][j-1];
				new_k += (k==0);
			    } else {
				cout << "found to edge " << it->second->name << " -> " << it->second->outNodes[k][j-1]->name << " " << k << " for read " << *id << endl;
				found_to = true;
				to = it->second->outNodes[k][j-1];
				new_k += 2*(k==3);
			    }
			}

			// either way, delete the edge 
			//delete_edge(it->second, it->second->outNodes[k][j-1], k);
		    }
		}
		if (found_from == true and found_to == true)
		{
		    cout << "add edge from  " << from->name << " -> " << to->name << " " << new_k << " for read  " << *id << endl;
		    add_edge(from->id, to->id, new_k, *id);
		}
	    }

	    // delete edges
	    for (uint k=0; k<4; ++k)
            {
                for (uint j=it->second->outNodes[k].size(); j!=0; --j)
                {
		    delete_edge(it->second, it->second->outNodes[k][j-1], k);
		}
	    }

	    // delete the node
	    cout << "delete node " << it->second->name << endl;
	    delete it->second;
            nodes.erase(it++);
	}
    }
}

void PanGraph::remove_low_covg_edges(const uint& thresh)
{
    cout << now() << "Edge cleaning with threshold " << thresh << endl;

    // remove all edges with less than thresh% coverage (mostly we want to remove all the singley occuring edges)
    for(map<uint32_t, PanNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        for (uint k=0; k<4; ++k)
        {
            for (uint32_t j=it->second->outNodes[k].size(); j!=0; --j)
            {
                if (it->second->outNodeCounts[k][it->second->outNodes[k][j-1]->id].size() < thresh)
                {
		    it->second->outNodes[k][j-1]->outNodeCounts[rev_orient(k)].erase(it->second->id);
                    it->second->outNodeCounts[k].erase(it->second->outNodes[k][j-1]->id);

		    cout << now() << "delete edge " << it->second->outNodes[k][j-1]->name << "->" << it->second->name << " " << rev_orient(k) << endl;
    		    it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].erase(std::remove(it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].begin(),
                                                  					    it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].end(),
                                                 					    it->second),
                                                                                it->second->outNodes[k][j-1]->outNodes[rev_orient(k)].end());
    		    cout << now() << "delete edge " << it->second->name << "->" << it->second->outNodes[k][j-1]->name << " " << k << endl;
                    it->second->outNodes[k].erase(std::remove(it->second->outNodes[k].begin(),
                                                  			it->second->outNodes[k].end(),
                                                  			it->second->outNodes[k][j-1]),
                                      			    it->second->outNodes[k].end());
                }
		assert(j-1 <= it->second->outNodes[k].size() || assert_msg("something has gone wrong when deleting edge " << it->second->name << "->" << it->second->outNodes[k][j-1]->name << " " << k << ". Did the edge occur twice in this read?"));
            }
        }
    }
}

void PanGraph::remove_isolated_nodes()
{
    cout << now() << "remove isolated nodes" << endl;
    bool found;
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
}

/*void PanGraph::delete_node(PanNode* n)
{
    // delete edges to node
    for (uint orientation=0; orientation<4; ++orientation)
    {
        for (uint32_t j=n->outNodes[k].size(); j!=0; --j)
        {
	    n->outNodes[k][j-1]->outNodes[r_orientation].erase(std::remove( n->outNodes[k][j-1]->outNodes[r_orientation].begin(),
                                                                            n->outNodes[k][j-1]->outNodes[r_orientation].end(),
                                                                            n),
                                                               n->outNodes[k][j-1]->outNodes[r_orientation].end());
	    n->outNodes[k][j-1]->outNodeCounts[r_orientation].erase(std::remove( n->outNodes[k][j-1]->outNodeCounts[r_orientation].begin(),
                                                                                 n->outNodes[k][j-1]->outNodeCounts[r_orientation].end(),
                                                                                 n->id),
                                                               n->outNodes[k][j-1]->outNodeCounts[r_orientation].end());
	}
    }

    // delete node
     cout << "delete node " << n->name;
    delete n;
    nodes.erase(n);
}*/

void PanGraph::clean(const uint32_t& covg)
{
    uint thresh = 0.05*covg;
    check_graph_symmetry();
    remove_low_covg_nodes(thresh);
    check_graph_symmetry();
    remove_low_covg_edges(thresh);
    check_graph_symmetry();
    remove_isolated_nodes();
    check_graph_symmetry();
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
	//cout << it->second->id << endl;
	
        handle << "S\t" << it->second->name << "\t*" << endl; //\tRC:i:" << it->second->foundHits.size() * k << endl;
        for (uint32_t j=0; j<it->second->outNodes[3].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[3][j]->name << "\t+\t0M\tRC:i:" << it->second->outNodeCounts[3][it->second->outNodes[3][j]->id].size() << endl;
        }
        for (uint32_t j=0; j<it->second->outNodes[2].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[2][j]->name << "\t-\t0M\tRC:i:" << it->second->outNodeCounts[2][it->second->outNodes[2][j]->id].size() << endl;
        }
	for (uint32_t j=0; j<it->second->outNodes[1].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[1][j]->name << "\t-\t0M\tRC:i:" << it->second->outNodeCounts[1][it->second->outNodes[1][j]->id].size() << endl;
        }
	for (uint32_t j=0; j<it->second->outNodes[0].size(); ++j)
        {
            handle << "L\t" << it->second->name << "\t+\t" << it->second->outNodes[0][j]->name << "\t-\t0M\tRC:i:" << it->second->outNodeCounts[0][it->second->outNodes[0][j]->id].size() << endl;
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
