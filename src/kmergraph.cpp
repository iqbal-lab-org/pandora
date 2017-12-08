#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <limits>
#include <cstdio>      /* NULL */
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <utility> // pair
#include "utils.h"
#include "kmernode.h"
#include "kmergraph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

KmerGraph::KmerGraph() {
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
KmerGraph::KmerGraph(const KmerGraph &other) {
    next_id = other.next_id;
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    thresh = other.thresh;
    KmerNodePtr n;

    // create deep copies of the nodes, minus the edges
    for (const auto &node : other.nodes) {
        n = make_shared<KmerNode>(*node.second);
        nodes[node.first] = n;
    }

    // now need to copy the edges
    for (auto c : other.nodes) {
        for (uint j = 0; j < c.second->outNodes.size(); ++j) {
            add_edge(nodes[c.first], nodes[c.second->outNodes[j]->id]);
        }
    }
}

// Assignment operator
KmerGraph &KmerGraph::operator=(const KmerGraph &other) {
    // check for self-assignment
    if (this == &other)
        return *this;

    // first we need to deallocate for any nodes already got!
    /*for (auto c: nodes) {
        delete c.second;
    }*/
    nodes.clear();

    // shallow copy no pointers
    next_id = other.next_id;
    num_reads = other.num_reads;
    shortest_path_length = other.shortest_path_length;
    k = other.k;
    p = other.p;
    thresh = other.thresh;
    KmerNodePtr n;

    // deep copy the vector of node pointers, excluding edges
    for (const auto &node : other.nodes) {
        n = make_shared<KmerNode>(*node.second);
        nodes[node.first] = n;
    }

    // now need to copy the edges
    for (auto c : other.nodes) {
        for (uint j = 0; j < c.second->outNodes.size(); ++j) {
            add_edge(nodes[c.first], nodes[c.second->outNodes[j]->id]);
        }
    }

    return *this;
}

KmerGraph::~KmerGraph() {
    clear();
}

void KmerGraph::clear() {
    /*for (auto c: nodes) {
        delete c.second;
    }*/
    nodes.clear();
    assert(nodes.empty());
    next_id = 0;
    num_reads = 0;
    shortest_path_length = 0;
    k = 0;
    p = 1;
    thresh = -25;
}

KmerNodePtr KmerGraph::add_node(const Path &p) {
    for (auto c : nodes) {
        if (c.second->path == p) {
            return c.second;
        }
    }

    // if we didn't find an existing node
    KmerNodePtr n (make_shared<KmerNode>(next_id, p));
    nodes[next_id] = n;
    assert(k == 0 or p.length() == 0 or p.length() == k);
    if (k == 0 and p.length() > 0) {
        k = p.length();
    }
    next_id++;
    if (next_id == reserved_size) {
        reserved_size *= 2;
        nodes.reserve(reserved_size);
    }
    return n;
}

KmerNodePtr KmerGraph::add_node_with_kh(const Path &p, const uint64_t &kh, const uint8_t &num) {
    KmerNodePtr n = add_node(p);
    n->khash = kh;
    n->num_AT = num;
    assert(n->khash < std::numeric_limits<uint64_t>::max());
    return n;
}

condition::condition(const Path &p) : q(p) {};

bool condition::operator()(const pair<uint32_t, KmerNodePtr> &kn) const { return kn.second->path == q; }

void KmerGraph::add_edge(const Path &from, const Path &to) {
    assert(from < to || assert_msg(from << " is not less than " << to));
    if (from == to) {
        return;
    }

    auto from_it = find_if(nodes.begin(), nodes.end(), condition(from));
    auto to_it = find_if(nodes.begin(), nodes.end(), condition(to));
    assert(from_it != nodes.end() && to_it != nodes.end());

    if (find(from_it->second->outNodes.begin(), from_it->second->outNodes.end(), to_it->second) ==
        from_it->second->outNodes.end()) {
        from_it->second->outNodes.push_back(to_it->second);
        to_it->second->inNodes.push_back(from_it->second);
        //cout << "added edge from " << (*from_it)->id << " to " << (*to_it)->id << endl;
    }
}

void KmerGraph::add_edge(KmerNodePtr from, KmerNodePtr to) {
    assert(from->path < to->path || assert_msg(from->id << " is not less than " << to->id));

    if (find(from->outNodes.begin(), from->outNodes.end(), to) == from->outNodes.end()) {
        from->outNodes.push_back(to);
        to->inNodes.push_back(from);
        //cout << "added edge from " << from->id << " to " << to->id << endl;
    }
}

void KmerGraph::check() {
    if (sorted_nodes.empty()) {
        sort_topologically();
    }

    // should not have any leaves, only nodes with degree 0 are start and end
    for (auto c = sorted_nodes.begin(); c != sorted_nodes.end(); ++c) {
        assert(!(*c)->inNodes.empty() or (*c) == sorted_nodes[0] ||
               assert_msg("node" << **c << " has inNodes size " << (*c)->inNodes.size()));
        assert(!(*c)->outNodes.empty() or (*c) == sorted_nodes.back() || assert_msg(
                "node" << **c << " has outNodes size " << (*c)->outNodes.size() << " and isn't equal to back node "
                       << *sorted_nodes.back()));
        for (auto d: (*c)->outNodes) {
            assert((*c)->path < d->path || assert_msg((*c)->path << " is not less than " << d->path));
            assert(find(c, sorted_nodes.end(), d) != sorted_nodes.end() ||
                   assert_msg(d->id << " does not occur later in sorted list than " << (*c)->id));
        }
    }
}

void KmerGraph::sort_topologically() {
    sorted_nodes.reserve(nodes.size());
    for (unordered_map<uint32_t, KmerNodePtr>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
        sorted_nodes.push_back(it->second);
    }
    sort(sorted_nodes.begin(), sorted_nodes.end(), pCompKmerNode());
}

void KmerGraph::get_next(const uint16_t kmer_id, const uint8_t thresh, unordered_set<uint16_t>& next_ids, vector<deque<KmerNodePtr>>& next_paths)
{
    // walk back in the graph until get hit or start node
    deque<KmerNodePtr> v = {nodes[kmer_id]};
    deque<deque<KmerNodePtr>> current_paths;
    current_paths.push_back(v);
    uint8_t num_shared_read;

    while(current_paths.size() > 0 and current_paths.size() < 1000)
    {
        v = current_paths.front();
        current_paths.pop_front();

        for (auto k : v.back()->outNodes)
        {
            v.push_back(k);
	    if (k->id == nodes[nodes.size()-1]->id or kmer_id == 0)
	    {
		next_paths.push_back(v);
                next_ids.insert(k->id);
	    } else if (k->covg[0] + k->covg[1] >= max((uint8_t)1,thresh))
            {
		// check if there are any reads which have a hit for both prev_id and kmer_id
                num_shared_read = 0;
                for (auto r : covgs)
                {
                    if (r[0][kmer_id]+r[1][kmer_id]+r[0][k->id]+r[1][k->id] >= 2)
                    {
                        num_shared_read += 1;
                        if (num_shared_read >= thresh)
			{
			    break;
			}
                    }
                }
                if (num_shared_read >= thresh)
                {
		    next_paths.push_back(v);
                    next_ids.insert(k->id);
                } else if (next_paths.size()>0 and v.size()<=3*next_paths[0].size()) {
		    current_paths.push_back(v);
		}		    
            } else {
                current_paths.push_back(v);
            }
            v.pop_back();
        }
    }
}

/*void KmerGraph::extend_paths_back(vector<deque<KmerNodePtr>>& paths_to_extend, const vector<deque<KmerNodePtr>>& path_extensions)
{
    if (path_extensions.size() == 0)
    {
	return;
    }

    deque<KmerNodePtr> a_path;
    vector<deque<KmerNodePtr>> extended_paths;
    extended_paths.reserve(500);

    for (auto it=paths_to_extend.begin(); it!=paths_to_extend.end();)
    {
        for (const auto d_ : path_extensions)
        {
            if(it->front() == d_.back())
	    {
            	a_path = *it;
            	a_path.insert(a_path.begin(), d_.begin(), --d_.end());
		extended_paths.push_back(a_path);
	    }
        }
        it++;
    }
    paths_to_extend = extended_paths;
}*/

void KmerGraph::extend_paths_forward(vector<deque<KmerNodePtr>>& paths_to_extend, const vector<deque<KmerNodePtr>>& path_extensions)
{
    if (path_extensions.size() == 0)
    {
        return;
    }

    deque<KmerNodePtr> a_path;
    vector<deque<KmerNodePtr>> extended_paths;
    extended_paths.reserve(500);
    bool extended;
    for (auto it=paths_to_extend.begin(); it!=paths_to_extend.end();)
    {
	extended = false;
	cout << "extend path ";
	for (auto n : *it)
	{
	   cout << n->id << " ";
	}
	cout << endl;
        for (const auto d_ : path_extensions)
        {
	    cout << "with extension ";
	    for (auto n : d_)
            {
           	cout << n->id << " ";
            }
            cout << endl;
            if(it->back() == d_.front())
	    {
                a_path = *it;
                a_path.insert(a_path.end(), ++d_.begin(), d_.end());
		extended_paths.push_back(a_path);
		extended = true;
	    }
        }
	if (extended == false)
	{
	    cout << "did not extend" << endl;
	    extended_paths.push_back(*it);
	}
        it++;
    }
    paths_to_extend = extended_paths;
}

/*void KmerGraph::find_compatible_paths(const uint16_t read_id, vector<deque<KmerNodePtr>>& paths)
{
    assert(covgs.size()>=read_id);
    assert(covgs[read_id].size() == 2);

    vector<uint16_t> prev(covgs[read_id][0].size(), 0);
    vector<uint16_t> next(covgs[read_id][0].size(), covgs[read_id][0].size()-1);
    vector<deque<KmerNodePtr>> v;
    v.reserve(500000);
    vector<vector<deque<KmerNodePtr>>> prev_paths(covgs[read_id][0].size(), v);
    vector<vector<deque<KmerNodePtr>>> next_paths(covgs[read_id][0].size(), v);
    uint8_t strand = 2;
    set<uint> hits_to_cover;

    cout << read_id << " ";
    // choose a strand - only one should have any hits
    for (uint st=0; st<2; ++st) {
        for (uint i = 0; i<covgs[read_id][st].size(); ++i) {
            if (covgs[read_id][st][i] > 0) {
                strand = st;
                break;
            }
        }
    }
    // and catch if no hits
    if (strand == 2)
    {
        return;
    }
    cout << strand << " ";

    // walk from each hit until cover all paths to the nearest hit
    for (uint i=0; i<covgs[read_id][strand].size(); ++i)
    {
        if (covgs[read_id][strand][i] > 0)
        {
            get_prev(read_id, strand, i, prev[i], prev_paths[i]);
            get_next(read_id, strand, i, next[i], next_paths[i]);
            hits_to_cover.insert(i);
        }
    }

    // collect this information into a set of paths to return
    vector<deque<KmerNodePtr>> paths_in_progress;
    uint j,i;
    cout << "cover with hits" << endl;
    while (hits_to_cover.size() > 0)
    {
	i = *hits_to_cover.begin();
	hits_to_cover.erase(i);

        paths_in_progress = prev_paths[i];
        j = prev[i];
	cout << i << "->" << j << "->";
        while(j > 0)
        {
            extend_paths_back(paths_in_progress, prev_paths[j]);
            hits_to_cover.erase(j);
            j = prev[j];
	    cout << j << "->";
        }
	cout << "done" << endl;
        extend_paths_forward(paths_in_progress, next_paths[i]);
        j = next[i];
	cout << i << "->" << j << "->";
        while(j < next.size()-1)
        {
            extend_paths_forward(paths_in_progress, next_paths[j]);
            hits_to_cover.erase(j);
            j = next[j];
	    cout << j << "->";
        }
	cout << "done" << endl;
        paths.insert(paths.end(),paths_in_progress.begin(), paths_in_progress.end());
    }
    cout << endl;
}*/

void KmerGraph::find_compatible_paths(const uint8_t thresh, vector<deque<KmerNodePtr>>& paths)
{
    unordered_set<uint16_t> v;
    vector<unordered_set<uint16_t>> next(nodes.size(), v);
    vector<deque<KmerNodePtr>> w;
    w.reserve(500);
    vector<vector<deque<KmerNodePtr>>> next_paths(nodes.size(), w);
    set<uint> hits_to_cover; // redefine this with custom sort based on sorted_nodes

    // walk from each hit until cover all paths to the nearest hit
    for (uint i=0; i<nodes.size(); ++i)
    {
        if (nodes[i]->covg[0] + nodes[i]->covg[1] >= thresh or i==0)
        {
	    cout << i << " as covg " << nodes[i]->covg[0] + nodes[i]->covg[1] << endl;
            get_next(i, thresh, next[i], next_paths[i]);
            hits_to_cover.insert(i);
        }
    }

    // collect this information into a set of paths to return
    vector<deque<KmerNodePtr>> path_extensions;
    uint i;

    cout << "cover with hits" << endl;
    for (auto p : next_paths[0])
    {
	if (hits_to_cover.find(p.back()->id)!=hits_to_cover.end())
	{
	    paths.push_back(p);
	}
    }
    while (hits_to_cover.size() > 0)
    {
        i = *hits_to_cover.begin();
        hits_to_cover.erase(i);
        cout << "hit " << i << "->";
	
	path_extensions.clear();
	for (auto j : next[i])
	{
	    cout << " " << j;
	    path_extensions.insert(path_extensions.end(), next_paths[j].begin(), next_paths[j].end());
	}
        cout << endl;

	extend_paths_forward(paths, path_extensions);
    }
}

/*void KmerGraph::find_all_compatible_paths(vector<deque<KmerNodePtr>>& all_paths, vector<vector<pair<uint16_t,uint16_t>>>& path_hits)
{
    cout << now() << "Finding all paths compatible with reads for node " << endl; 
    // adds all the compatible paths for all reads to paths
    // for each read, path_hits[read_id] gives a vector of <path_id, num_hits> pairs
    // for paths compatible with that read
    vector<deque<KmerNodePtr>> current_paths;
    vector<pair<uint16_t,uint16_t>> read_path_hits;
    current_paths.reserve(500000);
    read_path_hits.reserve(500000);
    all_paths.reserve(500000);
    path_hits.reserve(500000);
    uint path_i, num_hits;

    // collect the paths
    for (uint16_t read = 0; read!=covgs.size(); ++read)
    {
	cout << ".";
        path_hits.push_back(read_path_hits);
        find_compatible_paths(read, current_paths);
	cout << "update all_paths and path_hits for " << current_paths.size() << " paths" << endl;
        for (auto p : current_paths)
        {
            auto it = find (all_paths.begin(), all_paths.end(), p);
            if (it == all_paths.end())
            {
                all_paths.push_back(p);
                path_i = all_paths.size() -1;
            } else {
                path_i = distance(all_paths.begin(), it);
            }
            num_hits = 0;
            for (auto k : p)
            {
                num_hits += covgs[read][0][k->id] + covgs[read][1][k->id];
            }
	    cout << "path id " << path_i << " has length " << p.size() << " and num_hits " << num_hits << endl;
            path_hits[read].push_back(make_pair(path_i, num_hits));
        }
	cout << "done" << endl;
        current_paths.clear();
    }
    assert(path_hits.size() == covgs.size());
    cout << endl << now() << "Found " << all_paths.size() << " compatible paths" << endl;
}*/

void KmerGraph::find_all_compatible_paths(vector<deque<KmerNodePtr>>& all_paths, vector<vector<pair<uint16_t,uint16_t>>>& path_hits, const uint8_t thresh)
{
    // adds all the compatible paths for all reads to paths
    // for each read, path_hits[read_id] gives a vector of <path_id, num_hits> pairs
    // for paths compatible with that read
    vector<pair<uint16_t,uint16_t>> read_path_hits;
    read_path_hits.reserve(500);
    all_paths.reserve(500);
    path_hits.reserve(500);
    uint num_hits;

    // collect the paths
    find_compatible_paths(thresh, all_paths);
    cout << now() << "Found " << all_paths.size() << " compatible paths for node" << endl;

    // fill in hit data for reads
    for (uint16_t read = 0; read!=covgs.size(); ++read)
    {
        path_hits.push_back(read_path_hits);
        for (uint i=0; i<all_paths.size(); ++i)
        {
            num_hits = 0;
            for (auto k : all_paths[i])
            {
                num_hits += covgs[read][0][k->id] + covgs[read][1][k->id];
            }
	    if (num_hits > thresh)
	    {
                path_hits[read].push_back(make_pair(i, num_hits));
	    }
        }
        //cout << "found " << path_hits[read].size() << " compatible paths for read " << read << endl;
    }
    assert(path_hits.size() == covgs.size());
}

void KmerGraph::set_p(const float e_rate) {
    assert(k != 0);
    assert(0 < e_rate and e_rate < 1);
    p = 1 / exp(e_rate * k);
    //cout << "using p: " << p << endl;
}

float KmerGraph::prob(uint j) {
    assert(num_reads != 0);
    return prob(j, num_reads);
}

float KmerGraph::prob(uint j, uint num) {
    //prob of node j where j is node id (hence pos in nodes)
    assert(p != 1);
    assert(j < nodes.size());
    if (sorted_nodes.empty() and !nodes.empty()) {
        sort_topologically();
        check();
    }

    //cout << "prob of node " << j << " given " << num << " reads covering and covg " << nodes[j]->covg[0] << " , " << nodes[j]->covg[1];
    float ret;
    if (j == sorted_nodes[0]->id or j == sorted_nodes.back()->id) {
        ret = 0; // is really undefined
    } else if (nodes[j]->covg[0] + nodes[j]->covg[1] > num) {
        // under model assumptions this can't happen, but it inevitably will, so bodge
        ret = lognchoosek2(nodes[j]->covg[0] + nodes[j]->covg[1], nodes[j]->covg[0], nodes[j]->covg[1]) +
              (nodes[j]->covg[0] + nodes[j]->covg[1]) * log(p / 2);
        // note this may give disadvantage to repeat kmers
    } else {
        ret = lognchoosek2(num, nodes[j]->covg[0], nodes[j]->covg[1]) +
              (nodes[j]->covg[0] + nodes[j]->covg[1]) * log(p / 2) +
              (num - (nodes[j]->covg[0] + nodes[j]->covg[1])) * log(1 - p);
    }
    //cout << " is " << ret << endl;
    return ret;
}

float KmerGraph::find_max_path(vector<KmerNodePtr> &maxpath) {
    // finds a max likelihood path

    // need to catch if p not asserted...
    assert(p < 1 || assert_msg("p was not set in kmergraph"));
    assert(num_reads > 0 || assert_msg("num_reads was not set in kmergraph"));
    //p = 1/exp(e_rate*k);
    //cout << " with parameters n: " << num_reads << " and p: " << p << endl;
    //cout << "Kmer graph has " << nodes.size() << " nodes" << endl;

    // need to catch if thesh not set too...

    check();

    // create vectors to hold the intermediate values
    vector<float> M(nodes.size(), 0); // max log prob pf paths from pos i to end of graph
    vector<int> len(nodes.size(), 0); // length of max log path from pos i to end of graph
    vector<uint> prev(nodes.size(), nodes.size() - 1); // prev node along path
    float max_mean;
    int max_len;

    for (uint j = nodes.size() - 1; j != 0; --j) {
        max_mean = numeric_limits<float>::lowest();
        max_len = 0; // tie break with longest kmer path
        //cout << "node " << j-1 << " has " << sorted_nodes[j-1]->outNodes.size() << " outnodes" << endl;
        for (uint i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if ((sorted_nodes[j - 1]->outNodes[i]->id == sorted_nodes.back()->id and thresh > max_mean + 0.000001) or
                (M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] >
                 max_mean + 0.000001) or
                (max_mean - M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id] <=
                 0.000001 and len[sorted_nodes[j - 1]->outNodes[i]->id] > max_len)) {
                M[sorted_nodes[j - 1]->id] = prob(sorted_nodes[j - 1]->id) + M[sorted_nodes[j - 1]->outNodes[i]->id];
                len[sorted_nodes[j - 1]->id] = 1 + len[sorted_nodes[j - 1]->outNodes[i]->id];
                prev[sorted_nodes[j - 1]->id] = sorted_nodes[j - 1]->outNodes[i]->id;
                //cout << sorted_nodes[j-1]->id << " path: " << sorted_nodes[j-1]->path << " has prob: " << prob(j-1) << "  M: " << M[sorted_nodes[j - 1]->id] << " len: " << len[sorted_nodes[j - 1]->id] << " prev: " << prev[sorted_nodes[j - 1]->id] << endl;
                if (sorted_nodes[j - 1]->outNodes[i]->id != sorted_nodes.back()->id) {
                    max_mean = M[sorted_nodes[j - 1]->outNodes[i]->id] / len[sorted_nodes[j - 1]->outNodes[i]->id];
                    max_len = len[sorted_nodes[j - 1]->outNodes[i]->id];
                    //  cout << " and new max_mean: " << max_mean;
                } else {
                    max_mean = thresh;
                }
                //cout << endl;
            }
        }
        //cout << sorted_nodes[j-1]->id << " path: " << sorted_nodes[j-1]->path << "  M: " << M[sorted_nodes[j-1]->id] << " len: " << len[sorted_nodes[j-1]->id] << " prev: " << prev[sorted_nodes[j-1]->id] << endl;
    }
    // remove the final length added for the null start node
    len[0] -= 1;

    // extract path
    uint prev_node = prev[sorted_nodes[0]->id];
    while (prev_node < sorted_nodes.size() - 1) {
        //cout << prev_node << "->";
        maxpath.push_back(nodes[prev_node]);
        prev_node = prev[prev_node];
    }
    //cout << endl;
    //cout << "len[0]: " << len[0] << " maxpath.size(): " << maxpath.size() << " maxpath.back()->id: " << maxpath.back()->id << endl;

    assert(len[0] > 0 || assert_msg("found no path through kmer prg"));
    return M[0] / len[0];
}

vector<vector<KmerNodePtr>> KmerGraph::find_max_paths(uint num) {

    // save original coverges so can put back at the end
    vector<uint> original_covgs0, original_covgs1;
    for (uint i = 0; i != nodes.size(); ++i) {
        original_covgs0.push_back(nodes[i]->covg[0]);
        original_covgs1.push_back(nodes[i]->covg[1]);
    }

    // find num max paths
    //cout << "expected covg " << (uint)(p*num_reads/num) << endl;
    vector<vector<KmerNodePtr>> paths;
    vector<KmerNodePtr> maxpath;
    find_max_path(maxpath);
    //uint min_covg;
    paths.push_back(maxpath);

    while (paths.size() < num) {
        for (uint i = 0; i != maxpath.size(); ++i) {
            maxpath[i]->covg[0] -= min(maxpath[i]->covg[0], (uint) (p * num_reads / num));
            maxpath[i]->covg[1] -= min(maxpath[i]->covg[1], (uint) (p * num_reads / num));
        }
        maxpath.clear();
        find_max_path(maxpath);
        paths.push_back(maxpath);
    }

    // put covgs back
    for (uint i = 0; i != nodes.size(); ++i) {
        nodes[i]->covg[0] = original_covgs0[i];
        nodes[i]->covg[1] = original_covgs1[i];
    }

    return paths;
}

float KmerGraph::find_min_path(vector<KmerNodePtr> &maxpath) {
    // finds a paths with best minimum probability

    // need to catch if p not asserted...

    if (sorted_nodes.empty()) {
        sort_topologically();
        check();
    }

    // create vectors to hold the intermediate values
    vector<float> M(sorted_nodes.size(), 0); // min log prob of best path from pos i to end of graph
    vector<int> len(sorted_nodes.size(), 0); // length of min log path from pos i to end of graph
    vector<uint> prev(sorted_nodes.size(), sorted_nodes.size() - 1); // prev node along path
    float best_min;
    int best_len;

    for (uint j = sorted_nodes.size() - 1; j != 0; --j) {
        best_min = numeric_limits<float>::lowest();
        best_len = 0; // tie break with longest kmer path
        for (uint i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if ((sorted_nodes[j - 1]->outNodes[i]->id == sorted_nodes.size() - 1 and thresh > best_min + 0.000001) or
                (M[sorted_nodes[j - 1]->outNodes[i]->id] > best_min + 0.000001) or
                (best_min - M[sorted_nodes[j - 1]->outNodes[i]->id] <= 0.000001 and
                 len[sorted_nodes[j - 1]->outNodes[i]->id] > best_len)) {
                M[j - 1] = min(prob(j - 1), M[sorted_nodes[j - 1]->outNodes[i]->id]);
                len[j - 1] = 1 + len[sorted_nodes[j - 1]->outNodes[i]->id];
                prev[j - 1] = sorted_nodes[j - 1]->outNodes[i]->id;
                if (sorted_nodes[j - 1]->outNodes[i]->id != sorted_nodes.size() - 1) {
                    best_min = M[sorted_nodes[j - 1]->outNodes[i]->id];
                    best_len = len[sorted_nodes[j - 1]->outNodes[i]->id];
                } else {
                    best_min = thresh;
                }
            }
        }
    }

    // extract path
    uint prev_node = prev[0];
    while (prev_node < sorted_nodes.size() - 1) {
        maxpath.push_back(sorted_nodes[prev_node]);
        prev_node = prev[prev_node];
    }

    return M[0];
}

vector<vector<KmerNodePtr>> KmerGraph::get_random_paths(uint num_paths) {
    // find a random path through kmergraph picking ~uniformly from the outnodes at each point
    vector<vector<KmerNodePtr>> rpaths;
    vector<KmerNodePtr> rpath;
    uint i;

    time_t now;
    now = time(nullptr);
    srand((unsigned int) now);

    if (!nodes.empty()) {
        for (uint j = 0; j != num_paths; ++j) {
            i = rand() % nodes[0]->outNodes.size();
            rpath.push_back(nodes[0]->outNodes[i]);
            while (rpath.back() != nodes[nodes.size() - 1]) {
                if (rpath.back()->outNodes.size() == 1) {
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

float KmerGraph::prob_path(const vector<KmerNodePtr> &kpath) {
    float ret_p = 0;
    for (uint i = 0; i != kpath.size(); ++i) {
        ret_p += prob(kpath[i]->id);
    }
    uint len = kpath.size();
    if (kpath[0]->path.length() == 0) {
        len -= 1;
    }
    if (kpath.back()->path.length() == 0) {
        len -= 1;
    }
    if (len == 0) {
        len = 1;
    }
    return ret_p / len;
}

float KmerGraph::prob_paths(const vector<vector<KmerNodePtr>> &kpaths) {
    if (kpaths.empty()) {
        return 0; // is this the correct default?
    }

    // collect how much coverage we expect on each node from these paths
    vector<uint> path_node_covg(nodes.size(), 0);
    for (uint i = 0; i != kpaths.size(); ++i) {
        for (uint j = 0; j != kpaths[i].size(); ++j) {
            path_node_covg[kpaths[i][j]->id] += 1;
        }
    }

    // now calculate max likelihood assuming independent paths
    float ret_p = 0;
    uint len = 0;
    for (uint i = 0; i != path_node_covg.size(); ++i) {
        if (path_node_covg[i] > 0) {
            //cout << "prob of node " << nodes[i]->id << " which has path covg " << path_node_covg[i] << " and so we expect to see " << num_reads*path_node_covg[i]/kpaths.size() << " times IS " << prob(nodes[i]->id, num_reads*path_node_covg[i]/kpaths.size()) << endl;
            ret_p += prob(nodes[i]->id, num_reads * path_node_covg[i] / kpaths.size());
            if (nodes[i]->path.length() > 0) {
                len += 1;
            }
        }
    }

    if (len == 0) {
        len = 1;
    }

    //cout << "len " << len << endl;
    return ret_p / len;
}

void KmerGraph::save_covg_dist(const string &filepath) {

    ofstream handle;
    handle.open(filepath);

    for (auto c : nodes) {
        handle << c.second->covg[0] << "," << c.second->covg[1] << "," << (unsigned) c.second->num_AT << " ";
    }
    handle.close();
    return;
}

uint KmerGraph::min_path_length() {
    if (shortest_path_length > 0) {
        return shortest_path_length;
    }

    if (sorted_nodes.size() == 0) {
        sort_topologically();
        check();
    }

    vector<uint> len(sorted_nodes.size(), 0); // length of shortest path from node i to end of graph
    for (uint j = sorted_nodes.size() - 1; j != 0; --j) {
        for (uint i = 0; i != sorted_nodes[j - 1]->outNodes.size(); ++i) {
            if (len[sorted_nodes[j - 1]->outNodes[i]->id] + 1 > len[j - 1]) {
                len[j - 1] = len[sorted_nodes[j - 1]->outNodes[i]->id] + 1;
            }
        }
    }
    shortest_path_length = len[0];
    return len[0];
}

void KmerGraph::save(const string &filepath) {
    ofstream handle;
    handle.open(filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for (auto c : nodes) {
        handle << "S\t" << c.second->id << "\t" << c.second->path << "\tFC:i:" << c.second->covg[0] << "\t" << "\tRC:i:"
               << c.second->covg[1] << endl;//"\t" << (unsigned)nodes[i].second->num_AT << endl;
        for (uint32_t j = 0; j < c.second->outNodes.size(); ++j) {
            handle << "L\t" << c.second->id << "\t+\t" << c.second->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

void KmerGraph::load(const string &filepath) {
    string line;
    vector<string> split_line;
    stringstream ss;
    uint32_t id, covg, from, to;
    Path p;
    KmerNodePtr n;

    ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 4);
                id = stoi(split_line[1]);
                ss << split_line[2];
                ss >> p;
                ss.clear();
                //add_node(p);
                n = make_shared<KmerNode>(id, p);
                nodes[id] = n;
                next_id++;
                if (k == 0 and p.length() > 0) {
                    k = p.length();
                }
                assert(n->id == id);
                covg = stoi(split(split_line[3], "FC:i:")[0]);
                n->covg[0] = covg;
                covg = stoi(split(split_line[4], "RC:i:")[0]);
                n->covg[1] = covg;
                if (split_line.size() >= 6) {
                    n->num_AT = stoi(split_line[5]);
                }
            }
        }
        myfile.clear();
        myfile.seekg(0, myfile.beg);
        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
                if (split_line[2] == split_line[4]) {
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
}

bool KmerGraph::operator==(const KmerGraph &y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;
    }

    // false if have different nodes
    for (auto c : nodes) {
        // if node not equal to a node in y, then false
        auto found = find_if(y.nodes.begin(), y.nodes.end(), condition(c.second->path));
        if (found == y.nodes.end()) {
            return false;
        }

        // if the node is found but has different edges, then false
        if (c.second->outNodes.size() != found->second->outNodes.size()) { return false; }
        if (c.second->inNodes.size() != found->second->inNodes.size()) { return false; }
        for (uint32_t j = 0; j != c.second->outNodes.size(); ++j) {
            spointer_values_equal<KmerNode> eq = {c.second->outNodes[j]};
            if (find_if(found->second->outNodes.begin(), found->second->outNodes.end(), eq) ==
                found->second->outNodes.end()) { return false; }
        }

    }
    // otherwise is true
    return true;
}

bool pCompKmerNode::operator()(KmerNodePtr lhs, KmerNodePtr rhs) {
    return (lhs->path) < (rhs->path);
}

std::ostream &operator<<(std::ostream &out, KmerGraph const &data) {
    for (const auto c: data.nodes) {
        out << *(c.second);
    }
    return out;
}
