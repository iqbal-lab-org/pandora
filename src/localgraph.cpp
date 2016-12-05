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

void LocalGraph::add_node (const uint32_t& id, const string& seq, const Interval& pos, uint32_t nested_level)
{
    assert(seq.length() == pos.length); 
    map<uint32_t, LocalNode*>::iterator it=nodes.find(id);
    if(it==nodes.end())
    {
        LocalNode *n;
        n = new LocalNode(seq, pos, id, nested_level);
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

void LocalGraph::write_gfa (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0" << endl;
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

void LocalGraph::add_read_support_node (LocalNode* n)
{
    n->supported_by_reads = true;
    /*map<uint32_t, LocalNode*>::iterator it=nodes.find(id);
    if(it==nodes.end())
    {
	cout << "Error, node " << id << " is not in this prg!"<< endl;
    } else {
	it->second->supported_by_reads = true;
    }*/
    return; 
}

void LocalGraph::add_read_support_edge (LocalNode* f, LocalNode* t)
{
    if (find(f->supported_outNodes.begin(), f->supported_outNodes.end(), t)==f->supported_outNodes.end())
    {
        f->supported_outNodes.push_back(t);
        f->supported_out_deg+=1;
        t->supported_in_deg+=1;
    }
    /*map<uint32_t, LocalNode*>::iterator from_it=nodes.find(from);
    map<uint32_t, LocalNode*>::iterator to_it=nodes.find(to);
    assert((from_it!=nodes.end()) && (to_it!=nodes.end())); // NODE_MISSING_ERROR
    if((from_it!=nodes.end()) && (to_it!=nodes.end()))
    {
        LocalNode *f = (nodes.find(from)->second);
        LocalNode *t = (nodes.find(to)->second);
	//assert(f->supported_by_reads==true && t->supported_by_reads==true);
	// if it is not already there...
	if (find(f->supported_outNodes.begin(), f->supported_outNodes.end(), t)==f->supported_outNodes.end())
	{
            f->supported_outNodes.push_back(t);
	    f->supported_out_deg+=1;
	    t->supported_in_deg+=1;
	}
        //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    }*/
    return;
}

vector<deque<LocalNode*>> LocalGraph::node_step_forwards(vector<deque<LocalNode*>>& in_node_paths)
{
    vector<deque<LocalNode*>> out_node_paths;
    for(uint32_t i=0; i<in_node_paths.size(); ++i)
    {
	for(uint32_t j=0; j<in_node_paths[i].back()->outNodes.size(); ++j)
	{
	    out_node_paths.push_back(in_node_paths[i]);
	    out_node_paths.back().push_back(in_node_paths[i].back()->outNodes[j]);
	}
    }
    return out_node_paths;
}

vector<deque<LocalNode*>> LocalGraph::node_step_back(vector<deque<LocalNode*>>& in_node_paths)
{
    vector<deque<LocalNode*>> out_node_paths;
    for(uint32_t i=0; i<in_node_paths.size(); ++i)
    {
	for(map<uint32_t, LocalNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
	{
	    for (uint32_t j=0; j<it->second->outNodes.size(); ++j)
	    {
		if(it->second->outNodes[j]==in_node_paths[i].front())
		{
                    out_node_paths.push_back(in_node_paths[i]);
                    out_node_paths.back().push_front(it->second);
		}
            }
        }
    }
    return out_node_paths;
}

void LocalGraph::infer_read_supported_graph()
{
    // first work out what edges are supported
    // start by adding any edges from one supported node to another supported node
    for(map<uint32_t, LocalNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
	if (it->second->supported_by_reads==true)
	{
            for (uint32_t j=0; j<it->second->supported_outNodes.size(); ++j)
            {
		if (it->second->supported_outNodes[j]->supported_by_reads==true)
		{
		    add_read_support_edge(it->second, it->second->supported_outNodes[j]);
		}
            }
	}
    }

    // then look for nodes which have in degree or out degree < 1 and try to connect to nearest (in node steps) supported node
    for(map<uint32_t, LocalNode*>::iterator it=nodes.begin(); it!=nodes.end(); ++it)
    {
        if (it->second->supported_by_reads==true && it->second->supported_out_deg < 1)
	{
	    bool found_supported_node = false;
	    vector<deque<LocalNode*>> v = {{it->second}};
	    while(found_supported_node == false)
	    {
		v = node_step_forwards(v);
		for(uint32_t i=0; i<v.size(); ++i)
		{
		    if (v[i].back()->supported_by_reads==true or v[i].back()->id==nodes.size()-1) // ie if found a supported read, or reached the end of the graph
		    {
			// for each pair in node-path v[i] add it as a supported path
			for (uint32_t j = 1; j!= v[i].size(); ++j)
			{
			    add_read_support_node(v[i][j]);
			    add_read_support_edge(v[i][j-1], v[i][j]);
			}
			found_supported_node = true;
		    }
		}
	    }
	} else if (it->second->supported_by_reads==true && it->second->supported_in_deg < 1)
	{
	    bool found_supported_node = false;
            vector<deque<LocalNode*>> v = {{it->second}};
            while(found_supported_node == false)
            {
                v = node_step_back(v);
                for(uint32_t i=0; i<v.size(); ++i)
                {
                    if (v[i].front()->supported_by_reads==true or v[i].front()->id==0) // ie if found a supported read, or reached the front of the graph
                    {
                        // for each pair in node-path v[i] add it as a supported path
                        for (uint32_t j = 0; j!= v[i].size() - 1; ++j)
                        {
                            add_read_support_node(v[i][j]);
                            add_read_support_edge(v[i][j], v[i][j+1]);
                        }
                        found_supported_node = true;
                    }   
                }
            }
	}
    }
    return;
}

vector<deque<LocalNode*>> LocalGraph::get_read_supported_graph_paths(deque<LocalNode*>& in_node_path)
{
    vector<deque<LocalNode*>> out_node_paths;
    for(uint32_t j=0; j<in_node_path.back()->outNodes.size(); ++j)
    {
	if (in_node_path.back()->outNodes[j]->supported_by_reads==true)
	{
            out_node_paths.push_back(in_node_path);
            out_node_paths.back().push_back(in_node_path.back()->outNodes[j]);
	}
    }
    return out_node_paths;
}

void LocalGraph::write_read_supported_graph_paths(const string& filepath)
{
    ofstream handle;
    handle.open (filepath);

    uint32_t path_num = 0;
    vector<deque<LocalNode*>> v, new_supported_node_paths, supported_node_paths = {{nodes.begin()->second}};
    while (supported_node_paths.size()>0)
    {
	new_supported_node_paths.clear();
	for (uint32_t i=0; i!=supported_node_paths.size(); ++i)
	{
	    if(supported_node_paths[i].back()->id==nodes.size()-1)
	    {
		// write to file
                handle << ">path_" << path_num << endl;
		for (uint32_t j=0; j!=supported_node_paths[i].size(); ++j)
		{
		    handle << supported_node_paths[i][j]->seq;
		}
		handle << endl;
		path_num += 1;
	    } else {
		//step forward
		v = get_read_supported_graph_paths(supported_node_paths[i]);
	        new_supported_node_paths.insert(new_supported_node_paths.end(), v.begin(), v.end());
	    }
	}
    }
    return;
}
