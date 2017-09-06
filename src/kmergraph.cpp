#include <iostream>
#include <sstream>
#include <cstring>
#include <fstream>
#include <cassert>
#include <vector>
#include <limits>
#include <stdio.h>      /* NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
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
    thresh = -25;
}

// copy constructor
KmerGraph::KmerGraph(const KmerGraph& other)
{
    next_id = other.next_id;
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    thresh = other.thresh;

    // create deep copies of the nodes, minus the edges
    KmerNode* n;
    for (uint i=0; i<other.nodes.size(); ++i)
    {   
        assert(other.nodes[i]->id == i);
        n = new KmerNode(*other.nodes[i]);
        nodes.push_back(n);
    }

    // now need to copy the edges
    for (uint i=0; i<other.nodes.size(); ++i)
    {
	for (uint j=0; j<other.nodes[i]->outNodes.size(); ++j)
	{
	    add_edge(nodes[i], nodes[other.nodes[i]->outNodes[j]->id]);
	}
    }
}

// Assignment operator
KmerGraph& KmerGraph::operator=(const KmerGraph& other)
{
    // check for self-assignment
    if (this == &other)
        return *this;

    // first we need to deallocate for any nodes already got!
    for (auto c: nodes)
    {
        delete c;
    }
    nodes.clear();

    // shallow copy no pointers
    next_id = other.next_id;
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    thresh = other.thresh;

    // deep copy the vector of node pointers, excluding edges
    KmerNode* n;
    for (uint i=0; i<other.nodes.size(); ++i)
    {
	assert(other.nodes[i]->id == i);
	n = new KmerNode(*other.nodes[i]);
	nodes.push_back(n);
    }

    // now need to copy the edges
    for (uint i=0; i<other.nodes.size(); ++i)
    {
        for (uint j=0; j<other.nodes[i]->outNodes.size(); ++j)
        {
            add_edge(nodes[i], nodes[other.nodes[i]->outNodes[j]->id]);
        }
    }

    return *this;
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
    thresh = -25;
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

    if ( find(from->outNodes.begin(), from->outNodes.end(), to) == from->outNodes.end() )
    {
        from->outNodes.push_back(to);
        to->inNodes.push_back(from);
        //cout << "added edge from " << from->id << " to " << to->id << endl;
    }
    return;
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
    //cout << "reallocate ids" << endl;
    for (uint i=0; i!=nodes.size(); ++i)
    {
	nodes[i]->id = i;
    }
    return;
}

void KmerGraph::set_p(const float e_rate)
{
    p = 1/exp(e_rate*k);
}

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
    } else {
        ret = lognchoosek2(num_reads, nodes[j]->covg[0], nodes[j]->covg[1]) + (nodes[j]->covg[0]+nodes[j]->covg[1])*log(p/2) + 
		(num_reads-(nodes[j]->covg[0]+nodes[j]->covg[1]))*log(1-p);
    }
    return ret;
}

float KmerGraph::find_max_path(vector<KmerNode*>& maxpath)
{
    // finds a max likelihood path

    // need to catch if p not asserted...
    assert( p<1 || assert_msg("p was not set in kmergraph"));
    //p = 1/exp(e_rate*k);
    //cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    //cout << "Kmer graph has " << nodes.size() << " nodes" << endl;
    
    // need to catch if thesh not set too...

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    vector<int> len(nodes.size(), 0); // length of max log path from pos i to end of graph
    vector<uint> prev(nodes.size(), nodes.size()-1); // prev node along path
    float max_mean;
    int max_len;

    for (uint j=nodes.size()-1; j!=0; --j)
    {
        max_mean = numeric_limits<float>::lowest();
        max_len = 0; // tie break with longest kmer path
        for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
        {
            if ((nodes[j-1]->outNodes[i]->id == nodes.size()-1 and thresh > max_mean + 0.000001) or 
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
		    max_mean = thresh;
		}
		//cout << endl;
            }
        }
        //cout << j-1 << " path: " << nodes[j-1]->path << "  M: " << M[j-1] << " len: " << len[j-1] << " prev: " << prev[j-1] << endl;
    }
    // remove the final length added for the null start node
    len[0] -= 1;

    // extract path
    uint prev_node = prev[0];
    while (prev_node < nodes.size() - 1)
    {
        //cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }
    //cout << endl;
    //cout << "len[0]: " << len[0] << " maxpath.size(): " << maxpath.size() << " maxpath.back()->id: " << maxpath.back()->id << endl;

    return M[0]/len[0];
}

float KmerGraph::find_min_path(vector<KmerNode*>& maxpath)
{
    // finds a paths with best minimum probability

    // need to catch if p not asserted...

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // min log prob of best path from pos i to end of graph
    vector<int> len(nodes.size(), 0); // length of min log path from pos i to end of graph
    vector<uint> prev(nodes.size(), nodes.size()-1); // prev node along path
    float best_min;
    int best_len;

    for (uint j=nodes.size()-1; j!=0; --j)
    {
        best_min = numeric_limits<float>::lowest();
        best_len = 0; // tie break with longest kmer path
        for (uint i=0; i!=nodes[j-1]->outNodes.size(); ++i)
        {
            if ((nodes[j-1]->outNodes[i]->id == nodes.size()-1 and thresh > best_min + 0.000001) or 
                (M[nodes[j-1]->outNodes[i]->id] > best_min + 0.000001) or
                (best_min - M[nodes[j-1]->outNodes[i]->id] <= 0.000001 and len[nodes[j-1]->outNodes[i]->id] > best_len))
            {
                M[j-1] = min(prob(j-1), M[nodes[j-1]->outNodes[i]->id]);
                len[j-1] = 1 + len[nodes[j-1]->outNodes[i]->id];
                prev[j-1] = nodes[j-1]->outNodes[i]->id;
                if (nodes[j-1]->outNodes[i]->id != nodes.size()-1)
                {
                    best_min = M[nodes[j-1]->outNodes[i]->id];
                    best_len = len[nodes[j-1]->outNodes[i]->id];
                } else {
                    best_min = thresh;
                }
            }
        }
    }

    // extract path
    uint prev_node = prev[0];
    while (prev_node < nodes.size() - 1)
    {
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }

    return M[0];
}

vector<vector<KmerNode*>> KmerGraph::get_random_paths(uint num_paths)
{
    // find a random path through kmergraph picking ~uniformly from the outnodes at each point
    vector<vector<KmerNode*>> rpaths;
    vector<KmerNode*> rpath;
    uint i;

    time_t now;
    now = time(0);
    srand((unsigned int)now);

    if (nodes.size() > 0)
    {
	for (uint j=0; j!=num_paths; ++j)
	{
	    i = rand() % nodes[0]->outNodes.size();
            rpath.push_back(nodes[0]->outNodes[i]);
	    while (rpath.back() != nodes[nodes.size()-1])
	    {
	        if (rpath.back()->outNodes.size() == 1)
	        {
		    rpath.push_back(rpath.back()->outNodes[0]);
	        } else {
  	            i = rand() % rpath.back()->outNodes.size();
	            rpath.push_back(rpath.back()->outNodes[i]);
	        }
	    }
	    rpath.pop_back();
	    rpaths.push_back(rpath);
	    rpath.clear();
	}
    }
    return rpaths;
}

float KmerGraph::prob_path(const vector<KmerNode*>& kpath)
{
    float ret_p = 0;
    for (uint i=0; i!=kpath.size(); ++i)
    {
        ret_p += prob(kpath[i]->id);
    }
    uint len = kpath.size();
    if (kpath[0] == nodes[0])
    {
	len -= 1;
    }
    if (kpath.back() == nodes.back())
    {
        len -= 1;
    }
    if (len == 0)
    {
	len = 1;
    } 
    return ret_p/len;
}

void KmerGraph::save_covg_dist(const string& filepath)
{

    ofstream handle;
    handle.open(filepath);

    for (uint j=1; j!=nodes.size()-1; ++j)
    {
        handle << nodes[j]->covg[0] << "," << nodes[j]->covg[1] << "," << (unsigned)nodes[j]->num_AT << " ";
    }
    handle.close();
    return;
}

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
        handle << "S\t" << nodes[i]->id << "\t" << nodes[i]->path << "\tFC:i:" << nodes[i]->covg[0] << "\t" << "\tRC:i:" << nodes[i]->covg[1] << endl;//"\t" << (unsigned)nodes[i]->num_AT << endl;
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
		covg = stoi(split(split_line[3], "FC:i:")[0]);
		nodes.back()->covg[0] = covg;
		covg = stoi(split(split_line[4], "RC:i:")[0]);
		nodes.back()->covg[1] = covg;
		if (split_line.size() >= 6)
		{
		    nodes.back()->num_AT = stoi(split_line[5]);
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
