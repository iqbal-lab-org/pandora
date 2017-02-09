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

    for (uint m=0; m!=k; ++m)
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
                {prgs.push_back(s);}
        }
        cout << now() <<  "Number of LocalPRGs read: " << prgs.size() << endl;
        myfile.close();
    } else {
        cerr << "Unable to open PRG file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

void save_LocalPRG_kmer_paths(vector<LocalPRG*>& prgs, const string& prgfile)
{
    cout << now() << "Saving PRG minimizers" << endl;
    ofstream handle;
    handle.open (prgfile + ".mini");

    for (uint32_t i = 0; i != prgs.size(); ++i)
    {
        handle << i;
        for (uint j = 0; j!=prgs[i]->kmer_paths.size(); ++j)
        {
            handle << "\t" << prgs[i]->kmer_paths[j];
        }
        handle << endl;

    }
    handle.close();
    cout << now() << "Finished saving " << prgs.size() << " entries to file" << endl;
    return;
}

void load_LocalPRG_kmer_paths(vector<LocalPRG*>& prgs, const string& prgfile)
{
    cout << now() << "Loading PRG minimizers" << endl;
    uint32_t key;
    int c;
    Path p;

    ifstream myfile (prgfile + ".mini");
    if (myfile.is_open())
    {
        myfile >> key;
        while (myfile.good())
        {
            c = myfile.peek();
            if (c == '\n')
            {
		myfile.ignore(1,'\n');
                myfile >> key;
		assert(key < prgs.size());
            } else {
		myfile.ignore(1,'\t');
                myfile >> p;
                prgs[key]->kmer_paths.push_back(p);
            }
        }
    } else {
        cerr << "Unable to open PRG minimizer file " << prgfile << ".mini" << endl;
        exit(1);
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
	//cout << (*it)->kmer << endl;
        if (idx->minhash.find((*it)->kmer) != idx->minhash.end())
        {
	    for (vector<MiniRecord>::iterator it2=idx->minhash[(*it)->kmer].begin(); it2!=idx->minhash[(*it)->kmer].end(); ++it2)
            {
                //cout << (*it)->kmer << " : ";
	        hits->add_hit(s.id, *it, &(*it2));
		hit_count += 1;
            }
        }
    }
    cout << now() << "Found " << hit_count << " hits found for read " << name << " so size of MinimizerHits is now " << hits->hits.size() << endl;
    return;
}

void infer_localPRG_order_for_reads(const vector<LocalPRG*>& prgs, MinimizerHits* minimizer_hits, PanGraph* pangraph, const int max_diff, const uint32_t k)
{
    // this step infers the gene order for a read
    // orders hits from a set of minimizer hits, clusters them, removes noise hits, and adds the inferred gene order to the pangraph
    set<set<MinimizerHit*, pComp>,clusterComp> clusters_of_hits;

    if (minimizer_hits->hits.size() == 0) {return;}

    // First cluster hits matching same localPRG, not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    set<MinimizerHit*, pComp>::iterator mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHit*, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    float pn;
    for (set<MinimizerHit*, pComp>::iterator mh_current = ++minimizer_hits->hits.begin(); mh_current != minimizer_hits->hits.end(); ++mh_current)
    {
        //cout << "Hit: " << **mh_previous << endl;
        if((*mh_current)->read_id!=(*mh_previous)->read_id or (*mh_current)->prg_id!=(*mh_previous)->prg_id or (*mh_current)->strand!=(*mh_previous)->strand or (abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff)
        {
            pn = p_null(prgs, current_cluster, k);
	    //cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
            if (pn < 0.001)
            {
                //cout << "Found cluster of size: " << current_cluster.size() << " for prg " << (*mh_previous)->prg_id << endl;
                clusters_of_hits.insert(current_cluster);
	    }
            current_cluster.clear();
            current_cluster.insert(*mh_current);
        } else {
            current_cluster.insert(*mh_current);
        }
        mh_previous = mh_current;
    }
    pn = p_null(prgs, current_cluster, k);
    //cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
    if (pn < 0.001)
    {
        clusters_of_hits.insert(current_cluster);
        //cout << "Found final cluster of size: " << current_cluster.size() << " for prg " << (*mh_previous)->prg_id << endl;
    }
    //cout << "Found " << clusters_of_hits.size() << " clusters" << endl;

    // Next order clusters, remove contained ones, and add inferred order to pangraph    
    if (clusters_of_hits.size() == 0) { return;}
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id, *c_previous);
    //cout << "first cluster added " << (*(*c_previous).begin())->prg_id << endl; 
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) &&  ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) && ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start) ) // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, so force the next cluster to have at least a hit which is outside this region
        {
	    //cout << "added cluster with id " << (*(*c_current).begin())->prg_id << endl;
            pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
	    pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        //} else {
        //    //cout << "Contained cluster not added to order" << endl;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id)
	{
	    // if we just started looking at hits for a new read, add the first cluster
	    pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
            c_previous = c_current;
	} else {
	    //cout << "Did not add cluster. Criteria which may have failed:" << endl;
	    //cout << "read_ids equal: " << (*(*c_current).begin())->read_id << "==" << (*(*c_previous).begin())->read_id << ", " << ((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) << endl;
	    //cout << "prg_ids different: " << (*(*c_current).begin())->prg_id << "!=" << (*(*c_previous).begin())->prg_id << ", " << ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) << endl;
	    //cout << "not contained: " << (*--(*c_current).end())->read_interval.start << ">" << (*--(*c_previous).end())->read_interval.start << ", " << ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start) << endl;
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

void update_localPRGs_with_hits(PanGraph* pangraph, MinimizerHits* mh, const vector<LocalPRG*>& prgs, const uint32_t k, const float& e_rate)
{
    cout << now() << "Update coverages within genes" << endl;
    update_covgs_from_hits(prgs, mh);
    cout << now() << "Infer overlapped PRG paths for each PRG present" << endl;
    for(map<uint32_t, PanNode*>::iterator pnode=pangraph->nodes.begin(); pnode!=pangraph->nodes.end(); ++pnode)
    {
            cout << now() << "Looking at PRG " << pnode->second->id << endl;
            prgs[pnode->second->id]->infer_most_likely_prg_paths_for_corresponding_pannode(pnode->second, k, e_rate);
    }
}

/*void update_covgs_from_hits(vector<LocalPRG*>& prgs, MinimizerHits* mhs)
{
    for (uint32_t i=0; i!= prgs.size(); ++i)
    {
        prgs[i]->get_covgs(mhs);
    }
    return;
}*/
void update_covgs_from_hits(const vector<LocalPRG*>& prgs, MinimizerHits* mhs)
{
    // note that within mhs, hits which map to same prg should be together
    for (set<MinimizerHit*, pComp>::iterator mh = mhs->hits.begin(); mh != mhs->hits.end(); ++mh)
    {
        //cout << "prg_id: " << (*mh)->prg_id << endl;
        for (uint32_t i=0; i!= prgs.size(); ++i)
	{
	    if (prgs[i]->id == (*mh)->prg_id)
	    {
	        prgs[i]->update_covg_with_hit(*mh);
		break;
	    }
	}
    }
    return;
}

float p_null(const vector<LocalPRG*>& prgs, set<MinimizerHit*, pComp>& cluster_of_hits, uint32_t k)
{
    //cout << "finding p_null for cluster of size |x|=" << cluster_of_hits.size() << endl;
    float p = 0;

    // Assumes only one PRG in the vector has the id, (or works out p_null for the first occurance)
    assert(cluster_of_hits.size() > 0);
    bool id_found = false;

    for (uint32_t i=0; i!= prgs.size(); ++i)
    {
        if (prgs[i]->id == (*cluster_of_hits.begin())->prg_id)
	{
            //p = (1 - pow(1 - pow(0.25, k), cluster_of_hits.size()))*(1 - pow(1 - pow(0.25, k), prgs[i]->kmer_paths.size()));
            p = pow(1 - pow(1 - pow(0.25, k), prgs[i]->kmer_paths.size()), cluster_of_hits.size());
	    id_found = true;
            //cout << "p_null = " << p << ", using |y|=" << prgs[i]->kmer_paths.size() << " and |x|=" << cluster_of_hits.size() << " for prg " << prgs[i]->id << " and strand " << (*cluster_of_hits.begin())->strand << endl;
            return p;
	}
    }
    // if we have got here, then it couldn't be found
    cout << "did not find id" << endl;
    assert(id_found == true);
    return p;
}

