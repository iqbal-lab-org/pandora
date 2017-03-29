#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "utils.h"
#include "index.h"
#include "localPRG.h"
#include "seq.h"
#include "pangraph.h"
#include "pannode.h"
#include "minihits.h"
#include "minihit.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

string now()
{
    time_t now;
    string dt;
    
    now = time(0);
    dt = ctime(&now);
    return dt.substr(0,dt.length()-1) + " ";
}

vector<string> split(const string& query, const string& d)
{
    vector<string> v;
    string::size_type k = 0;
    string::size_type j = query.find(d, k);
    while (j!=string::npos) {
	if (j > k)
        {
	    v.push_back(query.substr(k, j-k));
	}
        k = j + d.size();
        j = query.find(d, k);
    }
    if (k < query.length())
    {
	v.push_back(query.substr(k));
    }	
    return v;
}

char complement(char n)
{
    switch(n)
    {
    case 'A':
    case 'a':
        return 'T';
    case 'T':
    case 't':
        return 'A';
    case 'G':
    case 'g':
        return 'C';
    case 'C':
    case 'c':
        return 'G';
    }
    //assert(false);
    return 'N';
}

string rev_complement(string s)
{
    transform(begin(s), end(s), begin(s), complement);
    reverse(s.begin(), s.end());
    return s;
}

// note this function is to be replaced due to considering log probs
float lognchoosek (uint32_t n, uint32_t k)
{
    assert(n >= k || assert_msg("Currently the model assumes that the most a given kmer (defined by position) can occur is once per read, i.e. an error somewhere else in the read cannot result in this kmer. If you are getting this message, then you have evidence of violation of this assumption. Either try using a bigger k, or come up with a better model"));
    float total = 0;

    for (uint m=n; m!=n-k; --m)
    {
	total += log(m);
    }

    for (uint m=1; m<k; ++m)
    {
	total -= log(m+1);
    }

    return total;
}

void read_prg_file(vector<LocalPRG*>& prgs, const string& filepath)
{
    cout << now() << "Loading PRGs from file " << filepath << endl;

    uint32_t id = 0;
    string name, read, line;
    LocalPRG *s;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        uint i = 0;
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!name.empty() && !read.empty())
                {
                    cout << now() << "Found PRG " << name << endl;
                    s = new LocalPRG(id, name, read);
                    
                    if (s!=nullptr)
                    {   
			prgs.push_back(s);
			id++;
		    } else {
			cerr << "Failed to make LocalPRG for " << name << endl;	
			exit(1);
		    }
                }
                name.clear();
                read.clear();
                if (!line.empty())
                {
                    name = line.substr(1);
                }
            }
            else
            {
                read += line;
            }
	    i++;
        }
        // and last entry
        if (!name.empty() && !read.empty())
        {
            cout << now() << "Found PRG " << name << endl;
            s = new LocalPRG(id, name, read);
            if (s!=nullptr)
            {prgs.push_back(s);
	    } else {
                cerr << "Failed to make LocalPRG for " << name << endl;
                exit(1);
            }
        }
        cout << now() <<  "Number of LocalPRGs read: " << prgs.size() << endl;
        myfile.close();
    } else {
        cerr << "Unable to open PRG file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

void load_PRG_kmergraphs(vector<LocalPRG*>& prgs, const string& prgfile)
{
    for (uint i=0; i!=prgs.size(); ++i)
    {
	prgs[i]->kmer_prg.load(prgfile + "." + to_string(i) + ".gfa");
    }
    return;
}

void add_read_hits(const uint32_t id, const string& name, const string& seq, MinimizerHits* hits, Index* idx, const uint32_t w, const uint32_t k)
{
    uint32_t hit_count = 0;
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and adds hits to a global MinimizerHits object
    Seq s(id, name, seq, w, k);
    for(set<Minimizer*, pMiniComp>::iterator it = s.sketch.begin(); it != s.sketch.end(); ++it)
    {
        if (idx->minhash.find((*it)->kmer) != idx->minhash.end())
        {
	    for (vector<MiniRecord>::iterator it2=idx->minhash[(*it)->kmer]->begin(); it2!=idx->minhash[(*it)->kmer]->end(); ++it2)
            {
	        hits->add_hit(s.id, *it, &(*it2));
		hit_count += 1;
            }
        }
    }
    hits->sort();
    cout << now() << "Found " << hit_count << " hits found for read " << name << " so size of MinimizerHits is now " << hits->hits.size() << endl;
    return;
}

void infer_localPRG_order_for_reads(const vector<LocalPRG*>& prgs, MinimizerHits* minimizer_hits, PanGraph* pangraph, const int max_diff, const uint32_t k)
{
    // this step infers the gene order for a read and adds this to the pangraph
    // by defining clusters of hits, keeping those which are not noise and
    // then adding the inferred gene ordering
    set<set<MinimizerHit*, pComp>,clusterComp> clusters_of_hits;

    if (minimizer_hits->hits.size() == 0) {return;}

    // First define clusters of hits matching same localPRG, not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    set<MinimizerHit*, pComp>::iterator mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHit*, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    float pn;
    for (set<MinimizerHit*, pComp>::iterator mh_current = ++minimizer_hits->hits.begin(); mh_current != minimizer_hits->hits.end(); ++mh_current)
    {
	cout << **mh_current << endl;
        if((*mh_current)->read_id!=(*mh_previous)->read_id or (*mh_current)->prg_id!=(*mh_previous)->prg_id or (*mh_current)->strand!=(*mh_previous)->strand or (abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff)
        {
	    cout << "change cluster because strand change==" << ((*mh_current)->strand!=(*mh_previous)->strand) << " or read position changed more than max_diff==" << ((abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff) << endl;
	    // keep clusters which have a low enough probability of occuring by chance
            pn = p_null(prgs, current_cluster, k);
            if (pn < 0.001)
            {
                clusters_of_hits.insert(current_cluster);
	    }
            current_cluster.clear();
        }
	current_cluster.insert(*mh_current);
        mh_previous = mh_current;
    }
    // keep final cluster if it has a low enough probability of occuring by chance
    pn = p_null(prgs, current_cluster, k);
    if (pn < 0.001)
    {
        clusters_of_hits.insert(current_cluster);
    }

    // Next order clusters, remove contained ones, and add inferred order to pangraph    
    if (clusters_of_hits.size() == 0) { return;}
    // to do this consider pairs of clusters in turn
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id, *c_previous);
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) &&  ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) && ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start) ) // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, so force the next cluster to have at least a hit which is outside this region
        {
            pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
	    pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id)
	{
	    // if we just started looking at hits for a new read, add the first cluster
	    pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
            c_previous = c_current;
	}
    }
    return;
}

void pangraph_from_read_file(const string& filepath, MinimizerHits* mh, PanGraph* pangraph, Index* idx, const vector<LocalPRG*>& prgs, const uint32_t w, const uint32_t k, const int max_diff)
{
    string name, read, line;
    uint32_t id = 0;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!read.empty()) // ok we'll allow reads with no name, removed
                {
		    cout << now() << "Found read " << name << endl;
		    cout << now() << "Add read hits" << endl;
                    add_read_hits(id, name, read, mh, idx, w, k);
		    id++;
                }
                name.clear();
                read.clear();
                if (!line.empty())
                {
                    name = line.substr(1);
                }
            }
            else
            {
                read += line;
            }
        }
        // and last entry
        if (!read.empty()) // allow reads with no name
        {
	    cout << now() << "Found read " << name << endl;
	    cout << now() << "Add read hits" << endl;
            add_read_hits(id, name, read, mh, idx, w, k);
        }
        //cout << "Number of reads found: " << id+1 << endl;
        cout << now() << "Infer gene orders and add to PanGraph" << endl;
        infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, k);
        myfile.close();
    } else {
        cerr << "Unable to open read file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

void update_localPRGs_with_hits(PanGraph* pangraph, const vector<LocalPRG*>& prgs) //, const uint32_t k, const float& e_rate, bool output_p_dist)
{
    for(map<uint32_t, PanNode*>::iterator pnode=pangraph->nodes.begin(); pnode!=pangraph->nodes.end(); ++pnode)
    {
        cout << now() << "Update coverages for PRG " << pnode->second->id << endl;
	for (set<MinimizerHit*, pComp_path>::iterator mh = pnode->second->foundHits.begin(); mh != pnode->second->foundHits.end(); ++mh)
	{
	    prgs[pnode->second->id]->update_covg_with_hit(*mh);
	}
	prgs[pnode->second->id]->kmer_prg.num_reads = pnode->second->foundReads.size();
	cout << now() << "Added " << prgs[pnode->second->id]->num_hits[1] << " hits in the forward direction and " << prgs[pnode->second->id]->num_hits[0] << " hits in the reverse" << endl;
    }
}

float p_null(const vector<LocalPRG*>& prgs, set<MinimizerHit*, pComp>& cluster_of_hits, uint32_t k)
{
    // Assumes only one PRG in the vector has the id, (or works out p_null for the first occurrence)
    assert(cluster_of_hits.size() > 0);
    assert((*cluster_of_hits.begin())->prg_id < prgs.size());

    uint32_t i = (*cluster_of_hits.begin())->prg_id;
    float p = pow(1 - pow(1 - pow(0.25, k), prgs[i]->kmer_prg.nodes.size()-2), cluster_of_hits.size());
    cout << "found cluster of size " << cluster_of_hits.size() << " against prg " << i << " with pnull " << p << endl;

    return p;
}

