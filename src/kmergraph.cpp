#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cassert>
#include <vector>
#include "utils.h"
#include "kmernode.h"
#include "kmergraph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

KmerGraph::KmerGraph()
{
    nodes.reserve(5000);
    next_id = 0;
    num_reads = 0;
    k = 0; // nb the kmer size is determined by the first non-null node added
    p = 1;
}

KmerGraph::~KmerGraph()
{
    clear();
}

void KmerGraph::clear()
{
    for (auto c: nodes)
    {
        delete c;
    }
    nodes.clear();
    assert(nodes.size() == 0);
    next_id = 0;
    num_reads = 0;
    k = 0;
    p = 1;
}

void KmerGraph::add_node (const Path& p)
{
    KmerNode *n;
    n = new KmerNode(next_id, p);
    pointer_values_equal<KmerNode> eq = { n };
    if ( find_if(nodes.begin(), nodes.end(), eq) == nodes.end() )
    {
	nodes.push_back(n);
	cout << "added node " << *n;
	assert(k==0 or p.length==0 or p.length==k);
	if (k == 0 and p.length > 0)
	{
	    k = p.length;
	}  
	next_id++;
    } else {
	delete n;
    }
    return;
}

void KmerGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    if (from == to)
    {
        return;
    }
    assert(from <= nodes.size() && to <= nodes.size());
    if ( find(nodes[from]->outNodes.begin(), nodes[from]->outNodes.end(), nodes[to]) == nodes[from]->outNodes.end() )
    {
        nodes[from]->outNodes.push_back(nodes[to]);
    }
    if ( find(nodes[to]->inNodes.begin(), nodes[to]->inNodes.end(), nodes[from]) == nodes[to]->inNodes.end() )
    {
        nodes[to]->inNodes.push_back(nodes[from]);
    }
    cout << "added edge from  " << from << " to " << to << endl;
    return;
}

condition::condition(const Path& p): q(p) {};
bool condition::operator()(const KmerNode* kn) const { return kn->path == q; }

void KmerGraph::add_edge (const Path& from, const Path& to)
{
    if (from == to)
    {
	return;
    }
    vector<KmerNode*>::iterator from_it = find_if(nodes.begin(), nodes.end(), condition(from));
    vector<KmerNode*>::iterator to_it = find_if(nodes.begin(), nodes.end(), condition(to));
    assert(from_it != nodes.end() && to_it != nodes.end());

    if ( find((*from_it)->outNodes.begin(), (*from_it)->outNodes.end(), (*to_it)) == (*from_it)->outNodes.end() )
    {
        (*from_it)->outNodes.push_back(*to_it);
	(*to_it)->inNodes.push_back((*from_it));
    }

    cout << "added edge from " << (*from_it)->id << " to " << (*to_it)->id << endl;
    return;
}

void KmerGraph::check (uint num_minikmers)
{
    // should have a node for every minikmer found, plus a dummy start and end
    assert(num_minikmers == 0 or nodes.size() == num_minikmers + 2 || assert_msg("nodes.size(): " << nodes.size() << " and num minikmers: " << num_minikmers));

    // should not have any leaves, only nodes with degree 0 are start and end
    for (auto c: nodes)
    {
	assert(c->inNodes.size() > 0 or c->id == 0 || assert_msg("node" << *c << " has inNodes size " << c->inNodes.size()));
	assert(c->outNodes.size() > 0 or c->id == nodes.size() - 1 || assert_msg("node" << *c << " has outNodes size " << c->outNodes.size()));
    }
    return;
}

/*vector<KmerNode*> KmerGraph::get_node_order()
{
    int num_bubble_starts = 0, num_bubble_ends = 0;
    vector<KmerNode*> return_order;
    return_order.reserve(next_id);
    vector<vector<KmerNode*>> nodes_by_level(10, return_order);
    assert(nodes_by_level.size() == 10);

    for (uint i=0; i!=nodes.size(); ++i)
    {
	if (nodes[i]->inNodes.size() > 1)
	{
	    num_bubble_ends += 1;
	    //cout << "num_bubble_ends is now " << num_bubble_ends << endl;
	}
	//cout << i << " " << num_bubble_starts - num_bubble_ends << endl;
	assert(num_bubble_starts - num_bubble_ends >= 0);
	nodes_by_level[num_bubble_starts - num_bubble_ends].push_back(nodes[i]);
	if (nodes[i]->outNodes.size() > 1)
        {
            num_bubble_starts += 1;
	    //cout << "num_bubble_starts is now " << num_bubble_starts << endl;
        }
    }
    
    //cout << "and output new order..." << endl;
    for (uint i=nodes_by_level.size(); i!=0; --i)
    {
	//cout << return_order.size() << " ";
        return_order.insert(return_order.end(), nodes_by_level[i-1].begin(), nodes_by_level[i-1].end());
    }
    return return_order;
}*/

float KmerGraph::prob(uint j, int dir)
{
    return lognchoosek(num_reads, nodes[j]->covg[dir]) + nodes[j]->covg[dir]*log(p) + (num_reads-nodes[j]->covg[dir])*log(1-p);
}

float KmerGraph::find_max_path(int dir, float e_rate, vector<KmerNode*>& maxpath)
{
    cout << now() << "Find kmer max path for direction " << dir;
    // update global p
    p = 1/exp(e_rate*k);
    cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    vector<int> len(nodes.size(), 1); // length of max log path from pos i to end of graph
    vector<uint> prev(nodes.size(), nodes.size()-1); // prev node along path
    float max_mean;

    M[nodes.size()-1] = prob(nodes.size()-1, dir);
    len[nodes.size()-1] = 1;
    //cout << nodes.size()-1 << "  M: " << M[nodes.size()-1] << " len: " << len[nodes.size()-1] << " prev: " << prev[nodes.size()-1] << endl;
    
    for (uint j=nodes.size()-1; j!=0; --j)
    {
	max_mean = numeric_limits<float>::lowest(); 
	for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
	{
	    if (M[nodes[j-1]->outNodes[i]->id]/len[nodes[j-1]->outNodes[i]->id] > max_mean)
	    {
		M[j-1] = prob(j-1, dir) + M[nodes[j-1]->outNodes[i]->id];
		len[j-1] = 1 + len[nodes[j-1]->outNodes[i]->id];
		prev[j-1] = nodes[j-1]->outNodes[i]->id;
		max_mean = M[nodes[j-1]->outNodes[i]->id]/len[nodes[j-1]->outNodes[i]->id];
	    }
	}
	//cout << j-1 << "  M: " << M[j-1] << " len: " << len[j-1] << " prev: " << prev[j-1] << endl;
    }

    // extract path
    uint prev_node = prev[0];
    while (prev_node < nodes.size() - 1)
    {
	//cout << prev_node << "->";
	maxpath.push_back(nodes[prev_node]);
	prev_node = prev[prev_node];
    }
    //cout << endl;

    return M[0]/len[0];
}

void KmerGraph::save (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for(uint i=0; i!=nodes.size(); ++i)
    {
        handle << "S\t" << nodes[i]->id << "\t" << nodes[i]->path << "\tRC:i:" << nodes[i]->covg[0] << "," << nodes[i]->covg[0] << endl;
        for (uint32_t j=0; j<nodes[i]->outNodes.size(); ++j)
        {
            handle << "L\t" << nodes[i]->id << "\t+\t" << nodes[i]->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

void KmerGraph::load (const string& filepath)
{
    string line;
    vector<string> split_line;
    stringstream ss;
    uint32_t id, covg, from, to;
    Path p;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
	    if (line[0] == 'S')
            {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = stoi(split_line[1]);
		ss << split_line[2];
		ss >> p;
		ss.clear();
                add_node(p);
		assert(nodes.back()->id == id);
		covg = stoi(split(split(split_line[3], "RC:i:")[0], ",")[0]);
		nodes.back()->covg[0] = covg;
		covg = stoi(split(split(split_line[3], "RC:i:")[0], ",")[1]);
                nodes.back()->covg[1] = covg;
            }
	}
        myfile.clear();
        myfile.seekg(0, myfile.beg);
        while ( getline (myfile,line).good() )
        {
            if (line[0] == 'L')
            {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                if (split_line[2] == split_line[4])
                {
                    from = stoi(split_line[1]);
                    to = stoi(split_line[3]);
                } else {
                    from = stoi(split_line[3]);
                    to = stoi(split_line[1]);
                }
                add_edge(from, to);
            }
        }
    } else {
        cerr << "Unable to open kmergraph file " << filepath << endl;
        exit(1);
    }
    return;		
}

bool KmerGraph::operator == (const KmerGraph& y) const
{
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;}

    // false if have different nodes
    for (uint i=0; i!=nodes.size(); ++i)
    {
        // if node not equal to a node in y, then false
        pointer_values_equal<KmerNode> eq = { nodes[i] };
	vector<KmerNode*>::const_iterator found = find_if(y.nodes.begin(), y.nodes.end(), eq);
        if ( found == y.nodes.end() )
	{
            return false;
	}

	// if the node is found but has different edges, then false
	if (nodes[i]->outNodes.size() != (*found)->outNodes.size()) {return false;}
	if (nodes[i]->inNodes.size() != (*found)->inNodes.size()) {return false;}
	for (uint32_t j=0; j!=nodes[i]->outNodes.size(); ++j)
        {
            pointer_values_equal<KmerNode> eq2 = { nodes[i]->outNodes[j] };
            if ( find_if((*found)->outNodes.begin(), (*found)->outNodes.end(), eq2) == (*found)->outNodes.end() )
            {return false;}
        }
	
    }
    // otherwise is true
    return true;
}

std::ostream& operator<< (std::ostream & out, KmerGraph const& data) {
    for (const auto c: data.nodes)
    {
        out << *c;
    }
    return out ;
}
