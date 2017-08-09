#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <cassert>
#include <algorithm>
//#include "FastxParser.hpp"
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

float lognchoosek2 (uint32_t n, uint32_t k1, uint32_t k2)
{
    assert(n >= k1+k2 || assert_msg("Currently the model assumes that the most a given kmer (defined by position) can occur is once per read, i.e. an error somewhere else in the read cannot result in this kmer. If you are getting this message, then you have evidence of violation of this assumption. Either try using a bigger k, or come up with a better model"));
    float total = 0;

    for (uint m=n; m!=n-k1-k2; --m)
    {
	total += log(m);
    }

    for (uint m=1; m<k1; ++m)
    {
	total -= log(m+1);
    }

    for (uint m=1; m<k2; ++m)
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

void load_PRG_kmergraphs(vector<LocalPRG*>& prgs, const uint& w, const uint& k, const string& prgfile)
{
    string prefix = "";
    size_t pos = prgfile.find_last_of("/");
    if (pos != std::string::npos)
    {
	prefix += prgfile.substr(0, pos);
	prefix += "/";
    }
    //cout << "prefix for kmerprgs dir is " << prefix << endl; 
    for (uint i=0; i!=prgs.size(); ++i)
    {
	prgs[i]->kmer_prg.load(prefix + "kmer_prgs/" + prgs[i]->name + ".k" + to_string(k) + ".w" + to_string(w) + ".gfa");
    }
    return;
}

//void add_read_hits(const uint32_t id, const string& name, const string& seq, MinimizerHits* hits, Index* idx, const uint32_t w, const uint32_t k)
void add_read_hits(Seq* s, MinimizerHits* hits, Index* idx)
{
    //cout << now() << "Search for hits for read " << s->name << " which has sketch size " << s->sketch.size() << " against index of size " << idx->minhash.size() << endl;
    uint32_t hit_count = 0;
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and adds hits to a global MinimizerHits object
    //Seq s(id, name, seq, w, k);
    for(set<Minimizer*, pMiniComp>::const_iterator it = s->sketch.begin(); it != s->sketch.end(); ++it)
    {
        if (idx->minhash.find((*it)->kmer) != idx->minhash.end())
        {
	    for (uint j=0; j!= idx->minhash[(*it)->kmer]->size(); ++j)
            {
	        hits->add_hit(s->id, *it, &(idx->minhash[(*it)->kmer]->operator[](j)));
		hit_count += 1;
            }
	//} else {
	//    cout << "did not find minimizer " << (*it)->kmer << " in index" << endl;
	}
    }
    //hits->sort();
    cout << now() << "Found " << hit_count << " hits found for read " << s->name << " so size of MinimizerHits is now " << hits->hits.size() + hits->uhits.size() << endl;
    return;
}

void infer_localPRG_order_for_reads(const vector<LocalPRG*>& prgs, MinimizerHits* minimizer_hits, PanGraph* pangraph, const int max_diff, const uint min_cluster_size)
{
    // this step infers the gene order for a read and adds this to the pangraph
    // by defining clusters of hits, keeping those which are not noise and
    // then adding the inferred gene ordering
    //cout << "sort" << endl;
    minimizer_hits->sort();
    //cout << "end sort" << endl;
    set<set<MinimizerHit*, pComp>,clusterComp> clusters_of_hits;

    if (minimizer_hits->hits.size() == 0) {return;}

    // First define clusters of hits matching same localPRG, not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    set<MinimizerHit*, pComp>::iterator mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHit*, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    for (set<MinimizerHit*, pComp>::iterator mh_current = ++minimizer_hits->hits.begin(); mh_current != minimizer_hits->hits.end(); ++mh_current)
    {
        if((*mh_current)->read_id!=(*mh_previous)->read_id or (*mh_current)->prg_id!=(*mh_previous)->prg_id or (*mh_current)->strand!=(*mh_previous)->strand or (abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff)
        {
	    // keep clusters which cover at least 10% of the shortest kmer path
            if (current_cluster.size() > max(prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()/10, min_cluster_size))
            {
                clusters_of_hits.insert(current_cluster);
	    /*} else {
		cout << "rejected hits" << endl;
		for (set<MinimizerHit*, pComp>::iterator p=current_cluster.begin(); p!=current_cluster.end(); ++p)
                {
                    cout << **p << endl;
                }*/
	    }
            current_cluster.clear();
        }
	current_cluster.insert(*mh_current);
        mh_previous = mh_current;
    }
    // keep final cluster if it has a low enough probability of occuring by chance
    if (current_cluster.size() > max(prgs[(*mh_previous)->prg_id]->kmer_prg.min_path_length()/20, min_cluster_size))
    {
        clusters_of_hits.insert(current_cluster);
    /*} else {
	cout << "rejected hits" << endl;
        for (set<MinimizerHit*, pComp>::iterator p=current_cluster.begin(); p!=current_cluster.end(); ++p)
        {
            cout << **p << endl;
        }*/
    }

    // Next order clusters, choose between those that overlap by too much
    cout << now() << "Found " << clusters_of_hits.size() << " clusters of hits " << endl;
    if (clusters_of_hits.size() == 0) { return;}
    // to do this consider pairs of clusters in turn
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    /*cout << "first cluster" << endl;
    for (set<MinimizerHit*, pComp>::iterator p=c_previous->begin(); p!=c_previous->end(); ++p)
    {
	cout << **p << endl;
    }*/
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
	/*cout << "current cluster" << endl;
        for (set<MinimizerHit*, pComp>::iterator p=c_current->begin(); p!=c_current->end(); ++p)
        {
            cout << **p << endl;
        }*/
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) && // if on same read and either
	  ((((*(*c_current).begin())->prg_id == (*(*c_previous).begin())->prg_id) && // same prg, different strand
	  ((*(*c_current).begin())->strand != (*(*c_previous).begin())->strand)) or // or cluster is contained
	  ((*--(*c_current).end())->read_interval.start <= (*--(*c_previous).end())->read_interval.start))) // i.e. not least one hit outside overlap 
	 // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, but could also impose no more than 2k hits in overlap
        {
	    if (c_previous->size() >= c_current->size())
	    {
		clusters_of_hits.erase(c_current);
		//cout << "erase current" << endl;
		c_current = c_previous;
	    } else {
		clusters_of_hits.erase(c_previous);
		//cout << "erase previous" << endl;
	    }
	}
	c_previous = c_current;
    }

    // Add inferred order to pangraph    
    cout << now() << "After removing contained clusters, have " << clusters_of_hits.size() << " clusters of hits " << endl;
    if (clusters_of_hits.size() == 0) { return;}
    // to do this consider pairs of clusters in turn
    c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, prgs[(*(*c_previous).begin())->prg_id]->name, (*(*c_previous).begin())->read_id, *c_previous);
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id)
        {
            pangraph->add_node((*(*c_current).begin())->prg_id, prgs[(*(*c_current).begin())->prg_id]->name, (*(*c_current).begin())->read_id, *c_current);
            pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id)
        {
            // if we just started looking at hits for a new read, add the first cluster
            pangraph->add_node((*(*c_current).begin())->prg_id, prgs[(*(*c_current).begin())->prg_id]->name, (*(*c_current).begin())->read_id, *c_current);
            c_previous = c_current;
        }
    }
    return;
}

/*void pangraph_from_read_file_new(const string& filepath, MinimizerHits* mh, PanGraph* pangraph, Index* idx, const vector<LocalPRG*>& prgs, const uint32_t w, const uint32_t k, const int max_diff)
{
    Seq *s;
    s = new Seq(0, "null", "", w, k);
    if (s==nullptr)
    {
        cerr << "Failed to create new Seq, something must be dying " << endl;
        exit (EXIT_FAILURE);
    }
    uint32_t id = 0;

    vector<string> files = {filepath};
    size_t nt = 8;
    size_t np = 1;
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, nt, np);
    parser.start();

    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser.getReadGroup();

    while (parser.refill(rg)) {
        // Here, rg will contain a chunk of reads
        // we can process.
        for (auto& read : rg) {
            //cout << "read name: " << read.name << ", seq: " << read.seq << '\n';
	    cout << now() << "Found read " << read.name << endl;
	    s->initialize(id, read.name, read.seq, w, k);
            cout << now() << "Add read hits" << endl;
            add_read_hits(s, mh, idx);
            id++;
        }
    }
    delete s;
    cout << now() << "Infer gene orders and add to PanGraph" << endl;
    infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff);
    return;
}*/

void pangraph_from_read_file(const string& filepath, MinimizerHits* mh, PanGraph* pangraph, Index* idx, const vector<LocalPRG*>& prgs, const uint32_t w, const uint32_t k, const int max_diff, const uint min_cluster_size, const uint genome_size)
{
    string name, read, line;
    uint64_t covg = 0;
    Seq *s;
    s = new Seq(0, "null", "", w, k);
    if (s==nullptr)
    {
        cerr << "Failed to create new Seq, something must be dying " << endl;
        exit (EXIT_FAILURE);
    }
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
                    s->initialize(id, name, read, w, k);
		    covg += s->seq.length();
		    cout << now() << "Add read hits" << endl;
		    add_read_hits(s, mh, idx);
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
            s->initialize(id, name, read, w, k);
	    covg += s->seq.length();
            cout << now() << "Add read hits" << endl;
            add_read_hits(s, mh, idx);
        }
	covg = covg/genome_size;
	cout << now() << "Estimated coverage: " << covg << endl;
        //cout << "Number of reads found: " << id+1 << endl;
        cout << now() << "Infer gene orders and add to PanGraph" << endl;
        infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, min_cluster_size);
	cout << now() << "Pangraph has " << pangraph->nodes.size() << " nodes" << endl;
        pangraph->clean(covg);
        cout << now() << "After cleaning, pangraph has " << pangraph->nodes.size() << " nodes" << endl;
	delete s;
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
        cout << now() << "Update coverages for PRG " << prgs[pnode->second->id]->name << endl;
	for (set<MinimizerHit*, pComp_path>::iterator mh = pnode->second->foundHits.begin(); mh != pnode->second->foundHits.end(); ++mh)
	{
	    prgs[pnode->second->id]->update_covg_with_hit(*mh);
	}
	prgs[pnode->second->id]->kmer_prg.num_reads = pnode->second->foundReads.size();
	cout << now() << "Added " << prgs[pnode->second->id]->num_hits[1] << " hits in the forward direction and " << prgs[pnode->second->id]->num_hits[0] << " hits in the reverse" << endl;
    }
}

/*float p_null(const vector<LocalPRG*>& prgs, set<MinimizerHit*, pComp>& cluster_of_hits, uint32_t k, const vector<uint32_t>& sketch_sizes)
{
    // Assumes only one PRG in the vector has the id, (or works out p_null for the first occurrence)
    assert(cluster_of_hits.size() > 0);
    assert((*cluster_of_hits.begin())->prg_id < prgs.size());

    uint32_t i = (*cluster_of_hits.begin())->prg_id;
    uint32_t sketch_size = sketch_sizes[(*cluster_of_hits.begin())->read_id];
    //float p = pow(1 - pow(1 - pow(0.25, k), prgs[i]->kmer_prg.nodes.size()-2), cluster_of_hits.size());
    float p_x = 1 - pow(1 - pow(0.25, k), prgs[i]->kmer_prg.nodes.size()-2);
    float p_y = 1 - pow(1 - pow(0.25, k), sketch_size);
    float J_null = p_x * p_y / (p_x + p_y - p_x * p_y);
    if (cluster_of_hits.size() > 2)
    {
	cout << "sketch size for read " << (*cluster_of_hits.begin())->read_id << " was " << sketch_size << " giving p_y " << p_y <<endl;
	cout << "Jnull is " << J_null << endl;
    }
    float p = 1;
    for (uint j=0; j!=cluster_of_hits.size(); ++j)
    {
	p -= nchoosek(sketch_size, j)*pow(J_null, j)*pow(1-J_null, sketch_size - j);
    }
    cout << "found cluster of size " << cluster_of_hits.size() << " against prg " << i << " with pnull " << p << endl;

    return p;
}*/

/*uint32_t nchoosek(const uint32_t n, const uint32_t k)
{
    uint32_t ret = 1;
    assert (n>=k);
    //assert (k>=0);

    if (n == k or k == 0)
    {
	ret = 1;
    } else {
        for (uint i=n; i>k; --i)
	{
	    ret *=ret;
	}
	for (uint i=1; i<=n-k; ++i)
        {
	    assert (i>0);
            ret /=ret;
        }
    }
    return ret;
}*/

/*float p_one(const vector<LocalPRG*>& prgs, set<MinimizerHit*, pComp>& cluster_of_hits, uint32_t k, float e_rate)
{
    // Assumes only one PRG in the vector has the id, (or works out p_null for the first occurrence)
    assert(cluster_of_hits.size() > 0);
    assert((*cluster_of_hits.begin())->prg_id < prgs.size());

    uint32_t i = (*cluster_of_hits.begin())->prg_id;
    float p = // put something to do with error rate, prob of getting this many correct

    return p;
}*/
