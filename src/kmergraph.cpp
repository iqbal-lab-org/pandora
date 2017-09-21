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
    reserved_size = 60000;
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
    for (unordered_map<uint32_t, KmerNode*>::const_iterator it=other.nodes.begin(); it!=other.nodes.end(); ++it)
    {
	n = new KmerNode((*it->second));
	nodes[it->first] = n;
    }

    // now need to copy the edges
    for (auto c : other.nodes)
    {
	for (uint j=0; j<c.second->outNodes.size(); ++j)
	{
	    add_edge(nodes[c.first], nodes[c.second->outNodes[j]->id]);
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
        delete c.second;
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
    for (unordered_map<uint32_t, KmerNode*>::const_iterator it=other.nodes.begin(); it!=other.nodes.end(); ++it)
    {
        n = new KmerNode(*(it->second));
        nodes[it->first] = n;
    }

    // now need to copy the edges
    for (auto c : other.nodes)
    {
        for (uint j=0; j<c.second->outNodes.size(); ++j)
        {
            add_edge(nodes[c.first], nodes[c.second->outNodes[j]->id]);
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
        delete c.second;
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
    for (auto c : nodes)
    {
	if (c.second->path == p)
	{
	    return c.second;
	}
    }

    // if we didn't find an existing node
    KmerNode *n;
    n = new KmerNode(next_id, p);
    nodes[next_id] = n;
    assert(k==0 or p.length()==0 or p.length()==k);
    if (k == 0 and p.length() > 0)
    {
        k = p.length();
    }
    next_id++;
    if (next_id == reserved_size)
    {
	reserved_size *= 2;
        nodes.reserve(reserved_size);
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
bool condition::operator()(const pair<uint32_t,KmerNode*>& kn) const { return kn.second->path == q; }

void KmerGraph::add_edge (const Path& from, const Path& to)
{
    assert(from < to ||assert_msg(from << " is not less than " << to) );
    if (from == to)
    {
	return;
    }

    unordered_map<uint32_t, KmerNode*>::iterator from_it = find_if(nodes.begin(), nodes.end(), condition(from));
    unordered_map<uint32_t, KmerNode*>::iterator to_it = find_if(nodes.begin(), nodes.end(), condition(to));
    assert(from_it != nodes.end() && to_it != nodes.end());

    if ( find(from_it->second->outNodes.begin(), from_it->second->outNodes.end(), to_it->second) == from_it->second->outNodes.end() )
    {
        from_it->second->outNodes.push_back(to_it->second);
	to_it->second->inNodes.push_back(from_it->second);
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

void KmerGraph::check ()
{
    if (sorted_nodes.size() == 0)
    {
	sort_topologically();
    }

    // should not have any leaves, only nodes with degree 0 are start and end
    for (vector<KmerNode*>::iterator c=sorted_nodes.begin(); c!=sorted_nodes.end(); ++c)
    {
	assert((*c)->inNodes.size() > 0 or (*c) == sorted_nodes[0] || assert_msg("node" << **c << " has inNodes size " << (*c)->inNodes.size()));
	assert((*c)->outNodes.size() > 0 or (*c) == sorted_nodes.back() || assert_msg("node" << **c << " has outNodes size " << (*c)->outNodes.size() << " and isn't equal to back node " << *sorted_nodes.back() ));
	for (auto d: (*c)->outNodes)
	{
	    assert((*c)->path < d->path || assert_msg((*c)->path << " is not less than " << d->path));
	    assert(find(c, sorted_nodes.end(), d) != sorted_nodes.end() || assert_msg(d->id << " does not occur later in sorted list than " << (*c)->id));
	}
    }
    return;
}

void KmerGraph::sort_topologically()
{
    sorted_nodes.reserve(nodes.size());
    for(unordered_map<uint32_t,KmerNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it) 
    {
  	sorted_nodes.push_back(it->second);
    }
    sort(sorted_nodes.begin(), sorted_nodes.end(), pCompKmerNode());
}

void KmerGraph::set_p(const float e_rate)
{
    p = 1/exp(e_rate*k);
}

float KmerGraph::prob(uint j)
{
    if (sorted_nodes.size() == 0 and nodes.size() > 0)
    {
	sort_topologically();
        check();
    }
	
    float ret;
    if (j==sorted_nodes[0]->id or j==sorted_nodes.back()->id)
    {    ret = 0; // is really undefined
    } else if (nodes[j]->covg[0]+nodes[j]->covg[1] > num_reads)
    {
	// under model assumptions this can't happen, but it inevitably will, so bodge
	ret = lognchoosek2(nodes[j]->covg[0]+nodes[j]->covg[1], nodes[j]->covg[0], nodes[j]->covg[1]) + 
		(nodes[j]->covg[0]+nodes[j]->covg[1])*log(p/2);
        // note this may give disadvantage to repeat kmers
    } else {
        ret = lognchoosek2(num_reads, nodes[j]->covg[0], nodes[j]->covg[1]) + 
		(nodes[j]->covg[0]+nodes[j]->covg[1])*log(p/2) + 
		(num_reads-(nodes[j]->covg[0]+nodes[j]->covg[1]))*log(1-p);
    }
    return ret;
}

float KmerGraph::find_max_path(vector<KmerNode*>& maxpath)
{
    // finds a max likelihood path

    // need to catch if p not asserted...
    assert( p<1 || assert_msg("p was not set in kmergraph"));
    assert( num_reads > 0 || assert_msg("num_reads was not set in kmergraph"));
    //p = 1/exp(e_rate*k);
    //cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    //cout << "Kmer graph has " << nodes.size() << " nodes" << endl;
    
    // need to catch if thesh not set too...

    check();

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
	//cout << "node " << j-1 << " has " << sorted_nodes[j-1]->outNodes.size() << " outnodes" << endl;
        for (uint i=0; i!=sorted_nodes[j-1]->outNodes.size(); ++i)
        {
            if ((sorted_nodes[j-1]->outNodes[i]->id == sorted_nodes.back()->id and thresh > max_mean + 0.000001) or 
		(M[sorted_nodes[j-1]->outNodes[i]->id]/len[sorted_nodes[j-1]->outNodes[i]->id] > max_mean + 0.000001) or
                (max_mean - M[sorted_nodes[j-1]->outNodes[i]->id]/len[sorted_nodes[j-1]->outNodes[i]->id] <= 0.000001 and len[sorted_nodes[j-1]->outNodes[i]->id] > max_len))
            {
                M[sorted_nodes[j-1]->id] = prob(sorted_nodes[j-1]->id) + M[sorted_nodes[j-1]->outNodes[i]->id];
                len[sorted_nodes[j-1]->id] = 1 + len[sorted_nodes[j-1]->outNodes[i]->id];
                prev[sorted_nodes[j-1]->id] = sorted_nodes[j-1]->outNodes[i]->id;
		//cout << j-1 << " path: " << sorted_nodes[j-1]->path << " has prob: " << prob(j-1) << "  M: " << M[j-1] << " len: " << len[j-1] << " prev: " << prev[j-1];
		if (sorted_nodes[j-1]->outNodes[i]->id != sorted_nodes.back()->id)
		{
                    max_mean = M[sorted_nodes[j-1]->outNodes[i]->id]/len[sorted_nodes[j-1]->outNodes[i]->id];
		    max_len = len[sorted_nodes[j-1]->outNodes[i]->id];
		  //  cout << " and new max_mean: " << max_mean;
		} else {
		    max_mean = thresh;
		}
		//cout << endl;
            }
        }
        //cout << j-1 << " path: " << sorted_nodes[j-1]->path << "  M: " << M[sorted_nodes[j-1]->id] << " len: " << len[sorted_nodes[j-1]->id] << " prev: " << prev[sorted_nodes[j-1]->id] << endl;
    }
    // remove the final length added for the null start node
    len[0] -= 1;

    // extract path
    uint prev_node = prev[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1)
    {
        //cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }
    //cout << endl;
    //cout << "len[0]: " << len[0] << " maxpath.size(): " << maxpath.size() << " maxpath.back()->id: " << maxpath.back()->id << endl;

    assert(len[0] > 0 || assert_msg("found no path through kmer prg"));
    return M[0]/len[0];
}

vector<vector<KmerNode*>> KmerGraph::find_max_paths(uint num)
{
    vector<vector<KmerNode*>> paths;
    vector<KmerNode*> maxpath;
    find_max_path(maxpath);
    //uint min_covg;
    paths.push_back(maxpath);

    while (paths.size() < num)
    {
	/*min_covg = 10000;
	for (uint i=0; i!=maxpath.size(); ++i)
	{
	    min_covg = min(min_covg, min(maxpath[i]->covg[0], maxpath[i]->covg[1]));	    
	}*/
	for (uint i=0; i!=maxpath.size(); ++i)
	{
	    maxpath[i]->covg[0] -= min(maxpath[i]->covg[0], (uint)p*num_reads/num);
	    maxpath[i]->covg[1] -= min(maxpath[i]->covg[1], (uint)p*num_reads/num);
	}
        maxpath.clear();
	find_max_path(maxpath);
	paths.push_back(maxpath);
    }
    return paths;	
}

float KmerGraph::find_min_path(vector<KmerNode*>& maxpath)
{
    // finds a paths with best minimum probability

    // need to catch if p not asserted...

    if (sorted_nodes.size() == 0)
    {
        sort_topologically();
        check();
    }

    // create vectors to hold the intermediate values
    vector<float> M(sorted_nodes.size(), 0); // min log prob of best path from pos i to end of graph
    vector<int> len(sorted_nodes.size(), 0); // length of min log path from pos i to end of graph
    vector<uint> prev(sorted_nodes.size(), sorted_nodes.size()-1); // prev node along path
    float best_min;
    int best_len;

    for (uint j=sorted_nodes.size()-1; j!=0; --j)
    {
        best_min = numeric_limits<float>::lowest();
        best_len = 0; // tie break with longest kmer path
        for (uint i=0; i!=sorted_nodes[j-1]->outNodes.size(); ++i)
        {
            if ((sorted_nodes[j-1]->outNodes[i]->id == sorted_nodes.size()-1 and thresh > best_min + 0.000001) or 
                (M[sorted_nodes[j-1]->outNodes[i]->id] > best_min + 0.000001) or
                (best_min - M[sorted_nodes[j-1]->outNodes[i]->id] <= 0.000001 and len[sorted_nodes[j-1]->outNodes[i]->id] > best_len))
            {
                M[j-1] = min(prob(j-1), M[sorted_nodes[j-1]->outNodes[i]->id]);
                len[j-1] = 1 + len[sorted_nodes[j-1]->outNodes[i]->id];
                prev[j-1] = sorted_nodes[j-1]->outNodes[i]->id;
                if (sorted_nodes[j-1]->outNodes[i]->id != sorted_nodes.size()-1)
                {
                    best_min = M[sorted_nodes[j-1]->outNodes[i]->id];
                    best_len = len[sorted_nodes[j-1]->outNodes[i]->id];
                } else {
                    best_min = thresh;
                }
            }
        }
    }

    // extract path
    uint prev_node = prev[0];
    while (prev_node < sorted_nodes.size() - 1)
    {
        maxpath.push_back(sorted_nodes[prev_node]);
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
    if (kpath[0]->path.length() == 0)
    {
	len -= 1;
    }
    if (kpath.back()->path.length() == 0)
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

    for (auto c : nodes)
    {
        handle << c.second->covg[0] << "," << c.second->covg[1] << "," << (unsigned)c.second->num_AT << " ";
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

    if (sorted_nodes.size() == 0)
    {
        sort_topologically();
        check();
    }

    vector<uint> len(sorted_nodes.size(), 0); // length of shortest path from node i to end of graph
    for (uint j=sorted_nodes.size()-1; j!=0; --j)
    {
        for (uint i=0; i!=sorted_nodes[j-1]->outNodes.size(); ++i)
        {
	    if (len[sorted_nodes[j-1]->outNodes[i]->id] + 1 > len[j-1])
	    {
		len[j-1] = len[sorted_nodes[j-1]->outNodes[i]->id] + 1;
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
    for(auto c : nodes)
    {
        handle << "S\t" << c.second->id << "\t" << c.second->path << "\tFC:i:" << c.second->covg[0] << "\t" << "\tRC:i:" << c.second->covg[1] << endl;//"\t" << (unsigned)nodes[i].second->num_AT << endl;
        for (uint32_t j=0; j<c.second->outNodes.size(); ++j)
        {
            handle << "L\t" << c.second->id << "\t+\t" << c.second->outNodes[j]->id << "\t+\t0M" << endl;
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
                n = new KmerNode(id, p);
		nodes[id] = n;
		next_id++;
		if (k == 0 and p.length() > 0)
                {
                    k = p.length();
                }
		assert(n->id == id);
		covg = stoi(split(split_line[3], "FC:i:")[0]);
		n->covg[0] = covg;
		covg = stoi(split(split_line[4], "RC:i:")[0]);
		n->covg[1] = covg;
		if (split_line.size() >= 6)
		{
		    n->num_AT = stoi(split_line[5]);
		}
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
    for (auto c : nodes)
    {
        // if node not equal to a node in y, then false
	unordered_map<uint32_t, KmerNode*>::const_iterator found = find_if(y.nodes.begin(), y.nodes.end(), condition(c.second->path));
        if ( found == y.nodes.end() )
	{
            return false;
	}

	// if the node is found but has different edges, then false
	if (c.second->outNodes.size() != found->second->outNodes.size()) {return false;}
	if (c.second->inNodes.size() != found->second->inNodes.size()) {return false;}
	for (uint32_t j=0; j!=c.second->outNodes.size(); ++j)
        {   pointer_values_equal<KmerNode> eq = { c.second->outNodes[j] };
            if ( find_if(found->second->outNodes.begin(), found->second->outNodes.end(), eq) == found->second->outNodes.end() )
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
        out << *(c.second);
    }
    return out ;
}
