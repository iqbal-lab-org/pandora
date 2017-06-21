#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cassert>
#include <vector>
#include <limits>
//#include <algorithm>
#include "utils.h"
#include "kmernode.h"
#include "kmergraph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

KmerGraph::KmerGraph()
{
    nodes.reserve(60000);
    next_id = 0;
    num_reads = 0;
    shortest_path_length = 0;    
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
    shortest_path_length = 0;
    k = 0;
    p = 1;
}

KmerNode* KmerGraph::add_node (const Path& p)
{
    KmerNode *n;
    n = new KmerNode(next_id, p);
    pointer_values_equal<KmerNode> eq = { n };
    vector<KmerNode*>::iterator it = find_if(nodes.begin(), nodes.end(), eq);
    if ( it == nodes.end() )
    {
	nodes.push_back(n);
	//cout << "added node " << *n;
	assert(k==0 or p.length()==0 or p.length()==k);
	if (k == 0 and p.length() > 0)
	{
	    k = p.length();
	}  
	next_id++;
    } else {
	//cout << "node " << *n << " was duplicate" << endl;
	delete n;
	n = *it;
    }

    return n;
}

KmerNode* KmerGraph::add_node_with_kh (const Path& p, const uint64_t& kh, const uint8_t& num)
{
    KmerNode *n = add_node(p);
    n->khash = kh;
    n->num_AT = num;
    assert(n->khash < std::numeric_limits<uint64_t>::max());
    return n;
}
    

/*void KmerGraph::add_edge (const uint32_t& from, const uint32_t& to)
{
    assert(from < to);
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
}*/

condition::condition(const Path& p): q(p) {};
bool condition::operator()(const KmerNode* kn) const { return kn->path == q; }

void KmerGraph::add_edge (const Path& from, const Path& to)
{
    assert(from < to ||assert_msg(from << " is not less than " << to) );
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
	//cout << "added edge from " << (*from_it)->id << " to " << (*to_it)->id << endl;
    }

    return;
}

void KmerGraph::add_edge (KmerNode* from, KmerNode* to)
{
    assert(from->path < to->path ||assert_msg(from->id << " is not less than " << to->id) );
    /*if (from->id > to->id)
    {
	cout << "switch id for nodes " << from->id << " and " << to->id << endl;
	vector<KmerNode*> new_nodes;
	new_nodes.insert(new_nodes.begin(), nodes.begin(), nodes.begin()+to->id);
	new_nodes.insert(new_nodes.end(), nodes.begin()+to->id+1, nodes.begin()+from->id+1);
	new_nodes.insert(new_nodes.end(),nodes.begin()+to->id, nodes.begin()+to->id+1);
	new_nodes.insert(new_nodes.end(), nodes.begin()+from->id+1, nodes.end());
	assert(new_nodes.size() == nodes.size());
	nodes = new_nodes;
	for (uint i=to->id; i<=from->id+1; ++i)
	{
	    cout << nodes[i]->id << "->" << i << ", ";
	    nodes[i]->id = i;
	}
	cout << endl;
    }*/

    if ( find(from->outNodes.begin(), from->outNodes.end(), to) == from->outNodes.end() )
    {
        from->outNodes.push_back(to);
        to->inNodes.push_back(from);
        //cout << "added edge from " << from->id << " to " << to->id << endl;
    }
    return;
}

/*void KmerGraph::copy_innodes (const Path& from, const Path& to)
{
    vector<KmerNode*>::iterator from_it = find_if(nodes.begin(), nodes.end(), condition(from));
    assert(from_it != nodes.end());

    for (uint i=0; i!=(*from_it)->inNodes.size(); ++i)
    {
	add_edge((*from_it)->inNodes[i]->path, to);
    }
    return;
}*/

set<Path> KmerGraph::get_innodes (const Path& from)
{
    set<Path> return_paths;
    vector<KmerNode*>::iterator from_it = find_if(nodes.begin(), nodes.end(), condition(from));
    assert(from_it != nodes.end());

    for (uint i=0; i!=(*from_it)->inNodes.size(); ++i)
    {
	return_paths.insert((*from_it)->inNodes[i]->path);
    }
    return return_paths;
}

void KmerGraph::check (uint num_minikmers)
{
    // should have a node for every minikmer found, plus a dummy start and end
    assert(num_minikmers == 0 or nodes.size() == num_minikmers || assert_msg("nodes.size(): " << nodes.size() << " and num minikmers: " << num_minikmers));

    // should not have any leaves, only nodes with degree 0 are start and end
    for (auto c: nodes)
    {
	assert(c->inNodes.size() > 0 or c->id == 0 || assert_msg("node" << *c << " has inNodes size " << c->inNodes.size()));
	assert(c->outNodes.size() > 0 or c->id == nodes.size() - 1 || assert_msg("node" << *c << " has outNodes size " << c->outNodes.size()));
	for (auto d: c->outNodes)
	{
	    assert(c->path < d->path || assert_msg(c->path << " is not less than " << d->path));
	    assert(c->id < d->id || assert_msg(c->id << " is not less than " << d->id));
	}
    }
    return;
}

void KmerGraph::sort_topologically()
{
    sort(nodes.begin(), nodes.end(), pCompKmerNode());
    /*vector<uint> num_seen_edges(nodes.size(), 0);
    vector<KmerNode*> found_order;
    deque<KmerNode*> to_add = {nodes[0]};
    KmerNode* kn;

    while (to_add.size() > 0)
    {
	kn = to_add.front();
	found_order.push_back(kn);
	//cout << kn->id << " ";
	to_add.pop_front();

	for (uint i=0; i!=kn->outNodes.size(); ++i)
	{
	    num_seen_edges[kn->outNodes[i]->id]+=1;
	    if (kn->outNodes[i]->inNodes.size() == num_seen_edges[kn->outNodes[i]->id])
	    {
		to_add.push_back(kn->outNodes[i]);
	    }
	}
    }
    //cout << endl << "found list size: " << found_order.size() << " as compared to " << nodes.size() << endl;
    assert(found_order.size() == nodes.size());
    nodes = found_order;
    */
    //cout << "reallocate ids" << endl;
    for (uint i=0; i!=nodes.size(); ++i)
    {
	nodes[i]->id = i;
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

float KmerGraph::prob(uint j)
{
    float ret;
    if (j==0 or j==nodes.size()-1)
    {    ret = 0; // is really undefined
    } else if (nodes[j]->covg[0]+nodes[j]->covg[1] > num_reads)
    {
	// under model assumptions this can't happen, but it inevitably will, so bodge
	ret = lognchoosek2(nodes[j]->covg[0]+nodes[j]->covg[1], nodes[j]->covg[0], nodes[j]->covg[1]) + (nodes[j]->covg[0]+nodes[j]->covg[1])*log(p/2);
        // note this may give disadvantage to repeat kmers
	//cout << "j: " << j << " special prob " << ret << " since has " << nodes[j]->covg[0] << " and " << nodes[j]->covg[1] << " hits" << endl;
    } else {
        ret = lognchoosek2(num_reads, nodes[j]->covg[0], nodes[j]->covg[1]) + (nodes[j]->covg[0]+nodes[j]->covg[1])*log(p/2) + 
		(num_reads-(nodes[j]->covg[0]+nodes[j]->covg[1]))*log(1-p);
        //cout << "j: " << j << " normal prob " << ret << " since has " << nodes[j]->covg[0] << " and " << nodes[j]->covg[1] << " hits" << endl;
	//cout << lognchoosek(num_reads, nodes[j]->covg[0]+nodes[j]->covg[1]) << " + " << (nodes[j]->covg[0]+nodes[j]->covg[1]) << "*" <<log(p) << " + ";
        //cout << (num_reads-(nodes[j]->covg[0]+nodes[j]->covg[1])) << "*" << log(1-p) << " + " << lognchoosek(nodes[j]->covg[0]+nodes[j]->covg[1], nodes[j]->covg[0]);
	//cout << " + " << (nodes[j]->covg[0]+nodes[j]->covg[1]) << "*" << log(0.5) << endl;
    }
    return ret;
}

/*float KmerGraph::find_max_path_forward(int dir, float e_rate, vector<KmerNode*>& maxpath)
{
    cout << now() << "Find forward kmer max path for direction " << dir;
    // update global p
    p = 1/exp(e_rate*k);
    cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    vector<int> len(nodes.size(), 1); // length of max log path from pos i to end of graph
    vector<uint> prev(nodes.size(), 0); // prev node along path
    float max_mean;
    int max_len;

    //M[0] = 0;
    //len[0] = 1;
    
    for (uint j=1; j!=nodes.size(); j++)
    {
	//cout << j << " ";
	max_mean = numeric_limits<float>::lowest(); 
	max_len = 0; // tie break with longest kmer path
	for (uint i=0; i!=nodes[j]->inNodes.size(); ++i)
	{
	    if ((M[nodes[j]->inNodes[i]->id]/len[nodes[j]->inNodes[i]->id] > max_mean) or
		(M[nodes[j]->inNodes[i]->id]/len[nodes[j]->inNodes[i]->id] == max_mean and len[nodes[j]->inNodes[i]->id] > max_len))
	    {
		M[j] = prob(j, dir) + M[nodes[j]->inNodes[i]->id];
		len[j] = 1 + len[nodes[j]->inNodes[i]->id];
		prev[j] = nodes[j]->inNodes[i]->id;
		max_mean = M[nodes[j]->inNodes[i]->id]/len[nodes[j]->inNodes[i]->id];
		max_len = len[nodes[j]->inNodes[i]->id];
	    }
	}
	//cout << j << "  M: " << M[j] << " len: " << len[j] << " prev: " << prev[j] << endl;
    }

    // extract path
    uint prev_node = prev[nodes.size()-1];
    while (prev_node > 0)
    {
	cout << prev_node << "->";
	maxpath.push_back(nodes[prev_node]);
	prev_node = prev[prev_node];
    }
    cout << endl;

    reverse(maxpath.begin(), maxpath.end());
    return M[nodes.size()-1]/len[nodes.size()-1];
}*/

//float KmerGraph::find_max_path_backward(int dir, float e_rate, vector<KmerNode*>& maxpath)
float KmerGraph::find_max_path(float e_rate, vector<KmerNode*>& maxpath)
{
    //cout << "running UPDTED VERSION" << endl;
    // update global p
    p = 1/exp(e_rate*k);
    //cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    //cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    vector<int> len(nodes.size(), 0); // length of max log path from pos i to end of graph
    vector<uint> prev(nodes.size(), nodes.size()-1); // prev node along path
    float max_mean;
    int max_len;

    //M[nodes.size()-1] = prob(nodes.size()-1);
    //len[nodes.size()-1] = 1;
    //cout << nodes.size()-1 << "  M: " << M[nodes.size()-1] << " len: " << len[nodes.size()-1] << " prev: " << prev[nodes.size()-1] << endl;

    for (uint j=nodes.size()-1; j!=0; --j)
    {
        max_mean = numeric_limits<float>::lowest();
        max_len = 0; // tie break with longest kmer path
        for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
        {
	    //cout << i << " ";

	    //if (M[nodes[j-1]->outNodes[i]->id] > log(0.005))
	    //{
                //cout << j-1 << " path: " << nodes[j-1]->path << " consider outnode: " << nodes[j-1]->outNodes[i]->id << " which has M: " << M[nodes[j-1]->outNodes[i]->id] << " len: " << len[nodes[j-1]->outNodes[i]->id] << " and the current max_mean: " << max_mean << endl;
	    //}

            if ((nodes[j-1]->outNodes[i]->id == nodes.size()-1 and -25 > max_mean + 0.000001) or 
		(M[nodes[j-1]->outNodes[i]->id]/len[nodes[j-1]->outNodes[i]->id] > max_mean + 0.000001) or
                (max_mean - M[nodes[j-1]->outNodes[i]->id]/len[nodes[j-1]->outNodes[i]->id] <= 0.000001 and len[nodes[j-1]->outNodes[i]->id] > max_len))
            {
                M[j-1] = prob(j-1) + M[nodes[j-1]->outNodes[i]->id];
                len[j-1] = 1 + len[nodes[j-1]->outNodes[i]->id];
                prev[j-1] = nodes[j-1]->outNodes[i]->id;
		//cout << j-1 << " path: " << nodes[j-1]->path << " has prob: " << prob(j-1) << "  M: " << M[j-1] << " len: " << len[j-1] << " prev: " << prev[j-1];
		if (nodes[j-1]->outNodes[i]->id != nodes.size()-1)
		{
                    max_mean = M[nodes[j-1]->outNodes[i]->id]/len[nodes[j-1]->outNodes[i]->id];
		    max_len = len[nodes[j-1]->outNodes[i]->id];
		  //  cout << " and new max_mean: " << max_mean;
		} else {
		    max_mean = log(0.005);
		}
		//cout << endl;
            }
        }
        //cout << j-1 << " path: " << nodes[j-1]->path << "  M: " << M[j-1] << " len: " << len[j-1] << " prev: " << prev[j-1] << endl;
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

void KmerGraph::save_covg_dist(const string& filepath)
{

    ofstream handle;
    handle.open(filepath);

    for (uint j=1; j!=nodes.size()-1; ++j)
    {
        handle << nodes[j]->covg[0] << "," << nodes[j]->covg[1] << "," << nodes[j]->num_AT << " ";
    }
    handle.close();
    return;
}

/*float KmerGraph::find_max_path_coverage(int dir, float e_rate, vector<KmerNode*>& maxpath)
{
    // maximise based on total coverage
    cout << now() << "Find coverage kmer max path for direction " << dir;
    // update global p
    p = 1/exp(e_rate*k);
    cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // create vectors to hold the intermediate values
    vector<uint> M(nodes.size(), 0); // max total covg from start to pos i
    vector<uint> prev(nodes.size(), nodes.size()-1); // prev node along path
    uint max_covg;

    for (uint j=nodes.size()-1; j!=0; --j)
    {
        max_covg = 0;
        for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
        {
            if (M[nodes[j-1]->outNodes[i]->id] >= max_covg)
            {
                M[j-1] = nodes[j-1]->covg[dir] + M[nodes[j-1]->outNodes[i]->id];
                prev[j-1] = nodes[j-1]->outNodes[i]->id;
                max_covg = M[nodes[j-1]->outNodes[i]->id];
		if (max_covg > 0)
		{
		    cout << j-1 << " " << nodes[j-1]->path << "  M: " << M[j-1] << " prev: " << prev[j-1] << endl;
		}
            }
        }
	//cout << j-1 << " " << *nodes[j-1] << "  M: " << M[j-1] << " prev: " << prev[j-1] << endl;
    }

    // extract path
    uint prev_node = prev[0];
    float ret_prob = 0;
    while (prev_node < nodes.size() - 1)
    {
        cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
	ret_prob += prob(prev_node, dir);
    }
    cout << endl;

    assert(maxpath.size() > 0);
    return ret_prob/maxpath.size();
}*/

uint KmerGraph::min_path_length()
{
    if (shortest_path_length > 0)
    {
	return shortest_path_length;
    }

    vector<uint> len(nodes.size(), 0); // length of shortest path from pos i to end of graph
    for (uint j=nodes.size()-1; j!=0; --j)
    {
        for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
        {
	    if (len[nodes[j-1]->outNodes[i]->id] + 1 > len[j-1])
	    {
		len[j-1] = len[nodes[j-1]->outNodes[i]->id] + 1;
	    }
	}
    }
    shortest_path_length = len[0];
    return len[0];
}

void KmerGraph::save (const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for(uint i=0; i!=nodes.size(); ++i)
    {
        handle << "S\t" << nodes[i]->id << "\t" << nodes[i]->path << "\tRC:i:" << nodes[i]->covg[1] << "\t" << (unsigned)nodes[i]->num_AT << endl;// << "," << nodes[i]->covg[0] << endl;
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
    KmerNode* n;

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
                //add_node(p);
                n = new KmerNode(next_id, p);
		nodes.push_back(n);
		next_id++;
		if (k == 0 and p.length() > 0)
                {
                    k = p.length();
                }
		assert(nodes.back()->id == id);
		covg = stoi(split(split_line[3], "RC:i:")[0]);
		nodes.back()->covg[0] = covg;
		if (split_line.size() >= 5)
		{
		    nodes.back()->num_AT = stoi(split_line[4]);
		}
		//covg = stoi(split(split(split_line[3], "RC:i:")[0], ",")[0]);
		//nodes.back()->covg[0] = covg;
		//covg = stoi(split(split(split_line[3], "RC:i:")[0], ",")[1]);
                //nodes.back()->covg[1] = covg;
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
                //add_edge(from, to);
                nodes[from]->outNodes.push_back(nodes[to]);
		nodes[to]->inNodes.push_back(nodes[from]);
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

bool pCompKmerNode::operator()(KmerNode* lhs, KmerNode* rhs) {
        return (lhs->path)<(rhs->path);
}

std::ostream& operator<< (std::ostream & out, KmerGraph const& data) {
    for (const auto c: data.nodes)
    {
        out << *c;
    }
    return out ;
}
