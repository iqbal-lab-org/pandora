#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <cassert>
#include "utils.h"
#include "index.h"
#include "localPRG.h"
#include "seq.h"
#include "pangraph.h"
#include "pannode.h"
#include "minihits.h"
#include "minihit.h"

void index_prg_file(vector<LocalPRG*>& prgs, const string& filepath, Index* idx, const uint32_t w, const uint32_t k)
{
    time_t now;
    uint32_t id = 0;
    string name, read, line, dt, sdt;
    LocalPRG *s;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
	//cout << "Opened prg file: " << filepath << endl;
        uint i = 0;
        while ( getline (myfile,line).good() )
        {
	    //cout << "reading line " << i << endl;
            if (line.empty() || line[0] == '>' )
            {
		//cout << "line empty or starts with >" << endl;
                if (!name.empty() && !read.empty())
                {
		    now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
                    cout << sdt << " Found PRG " << name << endl;
                    s = new LocalPRG(id, name, read);
		    s->minimizer_sketch(idx, w, k);
                    
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
		    //cout << "new name " << name << endl;
                }
            }
            else
            {
                read += line;
		//cout << "read starts " << read.substr(0,3) << endl; 
            }
	    i++;
        }
        // and last entry
        if (!name.empty() && !read.empty())
        {
            now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
            cout << sdt << " Found PRG " << name << endl;
            s = new LocalPRG(id, name, read);
            s->minimizer_sketch(idx, w, k);
            if (s!=nullptr)
                {prgs.push_back(s);}
        }
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt <<  " Number of LocalPRGs added to Index: " << prgs.size() << endl;
	cout << sdt << " Number of keys in Index: " << idx->minhash.size() << endl;
        myfile.close();
    } else {
        cerr << "Unable to open PRG file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

void add_read_hits(const uint32_t id, const string& name, const string& seq, MinimizerHits* hits, Index* idx, const uint32_t w, const uint32_t k)
{
    time_t now;
    string dt, sdt;
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
    now = time(0);
    dt = ctime(&now);
    sdt = dt.substr(0,dt.length()-1);
    cout << sdt << " Found " << hit_count << " hits found for read " << name << endl;
    return;
}

void infer_localPRG_order_for_reads(const vector<LocalPRG*>& prgs, MinimizerHits* minimizer_hits, PanGraph* pangraph, const int max_diff, const uint32_t cluster_thresh, const uint32_t k)
{
    // this step infers the gene order for a read
    // orders hits from a set of minimizer hits, clusters them, removes noise hits, and adds the inferred gene order to the pangraph
    set<set<MinimizerHit*, pComp>,clusterComp> clusters_of_hits;

    if (minimizer_hits->hits.size() == 0) {return;}

    // First cluster hits matching same localPRG, not more than max_diff read bases from the last hit (this last bit is to handle repeat genes). 
    // Keep clusters with more than cluster_thresh hits.
    set<MinimizerHit*, pComp>::iterator mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHit*, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    float pn;
    for (set<MinimizerHit*, pComp>::iterator mh_current = ++minimizer_hits->hits.begin(); mh_current != minimizer_hits->hits.end(); ++mh_current)
    {
        //cout << "Hit: " << **mh_previous << endl;
        if((*mh_current)->read_id!=(*mh_previous)->read_id or (*mh_current)->prg_id!=(*mh_previous)->prg_id or (*mh_current)->strand!=(*mh_previous)->strand or (abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff)
        {
            //if (current_cluster.size() > cluster_thresh)
            pn = p_null(prgs, current_cluster, k);
	    cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
            if (pn < 0.05)
            {
                //cout << "Found cluster of size: " << current_cluster.size() << " ending with " << **mh_previous << ". Next hit " << **mh_current << endl;
                clusters_of_hits.insert(current_cluster);
            } else {
		//cout << "Cluster not added had size " << current_cluster.size() << " and ended with " << **mh_previous << endl;
	    }
            current_cluster.clear();
            current_cluster.insert(*mh_current);
        } else {
            current_cluster.insert(*mh_current);
        }
        mh_previous = mh_current;
    }
    //if (current_cluster.size() > cluster_thresh)
    pn = p_null(prgs, current_cluster, k);
    cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
    if (pn < 0.05)
    {
        clusters_of_hits.insert(current_cluster);
        //cout << "Found final cluster of size: " << current_cluster.size() << " ending with " << **mh_previous << endl;
    } else {
        //cout << "Final cluster not added had size " << current_cluster.size() << " and ended with " << **mh_previous << endl;
    }

    // Next order clusters, remove contained ones, and add inferred order to pangraph
    if (clusters_of_hits.size() == 0) { return;}
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id, *c_previous);
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) &&  ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) && ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start + k - 1) ) // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, so force the next cluster to have at least a hit which is outside this region
        {
            pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
	    pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        //} else {
        //    cout << "Contained cluster not added to order" << endl;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id)
	{
	    // if we just started looking at hits for a new read, add the first cluster
	    pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id, *c_current);
            c_previous = c_current;
	}
    }
    return;
}

void pangraph_from_read_file(const string& filepath, PanGraph* pangraph, Index* idx, const vector<LocalPRG*>& prgs, const uint32_t w, const uint32_t k, const int max_diff, const uint32_t cluster_thresh)
{
    time_t now;
    string name, read, line, dt, sdt;
    uint32_t id = 0;
    MinimizerHits* mh;
    mh = new MinimizerHits();

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
	//cout << "Opened read file: " << filepath << endl;
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!read.empty()) // ok we'll allow reads with no name, removed
                {
    		    now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
		    cout << sdt << " Found read " << name << endl;
                    mh = new MinimizerHits();
		    now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
		    cout << sdt << " Add read hits" << endl;
                    add_read_hits(id, name, read, mh, idx, w, k);
		    /*now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
		    cout << sdt << " Infer gene orders and add to PanGraph" << endl;
		    infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, cluster_thresh, k);
		    now = time(0);
                    dt = ctime(&now);
                    sdt = dt.substr(0,dt.length()-1);
		    cout << sdt << " Update coverages within genes" << endl;
		    update_covgs_from_hits(prgs, mh);
		    delete mh;*/
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
	    now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
	    cout << sdt << " Found read " << name << endl;
            mh = new MinimizerHits();
            now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
	    cout << sdt << " Add read hits" << endl;
            add_read_hits(id, name, read, mh, idx, w, k);
	    /*now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
	    cout << sdt << " Infer gene orders and add to PanGraph" << endl;
            infer_localPRG_order_for_read(prgs, mh, pangraph, max_diff, cluster_thresh, k);
            now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
	    cout << sdt << " Update coverages within genes" << endl;
	    update_covgs_from_hits(prgs, mh);
	    delete mh;*/
        }
        //cout << "Number of reads found: " << id+1 << endl;
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt << " Infer gene orders and add to PanGraph" << endl;
        infer_localPRG_order_for_reads(prgs, mh, pangraph, max_diff, cluster_thresh, k);
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt << " Update coverages within genes" << endl;
        update_covgs_from_hits(prgs, mh);
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt << " Infer overlapped PRG paths for each PRG present" << endl;
        for(map<uint32_t, PanNode*>::iterator pnode=pangraph->nodes.begin(); pnode!=pangraph->nodes.end(); ++pnode) 
	{
	    now = time(0);
            dt = ctime(&now);
            sdt = dt.substr(0,dt.length()-1);
            cout << sdt << " Looking at PRG " << pnode->second->id << endl;
	    infer_most_likely_prg_path_for_pannode(prgs, pnode->second, k, 0.0015);
	}
        delete mh;
        myfile.close();
    } else {
        cerr << "Unable to open read file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
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
    //uint32_t current_id = (*mhs->hits.begin())->prg_id;
    //deque<MinimizerHit*> current_hits;
    //cout << "prgs vector length: " << prgs.size() << endl;
    for (set<MinimizerHit*, pComp>::iterator mh = mhs->hits.begin(); mh != mhs->hits.end(); ++mh)
    {
        //cout << "prg_id: " << (*mh)->prg_id << endl;
	prgs[(*mh)->prg_id]->update_covg_with_hit(*mh);
	/*if ((*mh)->prg_id == current_id)
	{
	    current_hits.push_back(*mh);
	} else {
	    cout << "Update covg on " << current_id << " with " << current_hits.size() << " hits" << endl;
	    prgs[current_id]->update_covg_with_hits(current_hits);
	    current_hits.clear();
	    current_hits.push_back(*mh);
	    current_id = (*mh)->prg_id;
	}*/
    }
    // and last one
    //cout << "Update covg on " << current_id << " with " << current_hits.size() << " hits" << endl;
    //prgs[current_id]->update_covg_with_hits(current_hits);
    return;
}

float p_null(const vector<LocalPRG*>& prgs, set<MinimizerHit*, pComp>& cluster_of_hits, uint32_t k)
{
   float p = (1 - pow(1 - pow(0.25, k), cluster_of_hits.size()))*(1 - pow(1 - pow(0.25, k), prgs[(*cluster_of_hits.begin())->prg_id]->num_minis));
   return p;
}

void infer_most_likely_prg_path_for_pannode(const vector<LocalPRG*>& prgs, PanNode* pnode, uint32_t k, float e_rate)
{
    // note that within the foundHits, hits which map to same read, prg and strand will be grouped together
    // we want to know for each prg minimizer, how much support there is from the reads
    // start by counting how many hits against each minimizing_kmer of prg
    cout << "start by counting how many hits against each minimizing kmer in prg" << endl;
    vector<uint32_t> kmer_path_hit_counts(prgs[pnode->id]->kmer_paths.size(),0);  
    cout << "there are " << prgs[pnode->id]->kmer_paths.size() << " minimizing kmers to count hits against" << endl;
    set<MinimizerHit*, pComp_path>::iterator mh_previous = pnode->foundHits.begin();
    uint32_t num_hits = 1, num_kmers=0;
    for (set<MinimizerHit*, pComp_path>::iterator mh = (pnode->foundHits.begin()); mh != pnode->foundHits.end();)
    {
	mh++;
	if (mh==pnode->foundHits.end() or (*mh)->strand != (*mh_previous)->strand or !((*mh)->prg_path == (*mh_previous)->prg_path))
	{
	    for (uint32_t i=num_kmers; i!=prgs[pnode->id]->kmer_paths.size(); ++i)
	    {
		if (prgs[pnode->id]->kmer_paths[i]==(*mh_previous)->prg_path)
		{
		    //cout << "found path " << prgs[pnode->id]->kmer_paths[i] << " == " << (*mh_previous)->prg_path << endl;
		    kmer_path_hit_counts[i] = num_hits;
		    //cout << "kmer_path_hit_counts[" << i << "] = " << num_hits << endl;
		    num_hits = 1;
		    num_kmers +=1;
		    break;
		}
	    }
	    assert(num_hits==1);
	} else {
	    num_hits += 1;
	    //cout << num_hits << " ";
	}
	mh_previous++;
    }

    uint32_t sum_of_elems = 0;
    for (uint32_t n : kmer_path_hit_counts)
    {sum_of_elems += n;}
    cout << "after adding counts of the " << pnode->foundHits.size() << " hits, have still got " << sum_of_elems << " hits added to tallies" << endl;
    assert(sum_of_elems==pnode->foundHits.size());

    // now for each of the minimizing kmers, work out the prob of seeing this number of hits given the number of reads
    // this is the bit where I assume that we have an independent trial for each read (binomial hit counts for true kmers)
    cout << "next work out prob of seeing this number of hits against a kmer assuming it is truly present" << endl;
    float p_thresh=0.5;
    vector<float> kmer_path_probs;
    float p_kmer, p=1/exp(e_rate*k), p_max=0, p_min=1;
    uint32_t n = pnode->foundReads.size(), big_p_count=0;
    cout << "n: " << n << ", p: " << p << endl;
    cout << "count:0 " << nchoosek(n,0) << " * " << pow(p,0) << " * " << pow(1-p,n-0) << endl;
    cout << "count:1 " << nchoosek(n,1) << " * " << pow(p,1) << " * " << pow(1-p,n-1) << endl;
    for (uint32_t i=0; i!=kmer_path_hit_counts.size(); ++i)
    {
        p_kmer = nchoosek(n, kmer_path_hit_counts[i])*pow(p,kmer_path_hit_counts[i])*pow(1-p,n-kmer_path_hit_counts[i]);
        kmer_path_probs.push_back(p_kmer);
        p_max = max(p_max, p_kmer);
        p_min = min(p_min, p_kmer);
	if(p_kmer>p_thresh)
        {
            //cout << p_kmer << " ";
            big_p_count+=1;
        }
    }
    cout << endl;
    cout << "found " << big_p_count << " probs > " << p_thresh << " with max and min probs " << p_max << ", " << p_min << endl;

    // now work out which nodes of graph have support from a kmer with prob > p_thresh
    cout << "now look at which nodes of prg have support" << endl;
    uint32_t num_supported_nodes = 0;
    bool found_support;
    // i.e. for each prg node...
    for (map<uint32_t,LocalNode*>::iterator n=prgs[pnode->id]->prg.nodes.begin(); n!=prgs[pnode->id]->prg.nodes.end(); ++n)
    {
	found_support = false;
	uint32_t i=0;
	//...check each minimizer...
	while(found_support == false and i!=prgs[pnode->id]->kmer_paths.size())
	{
	    //...and if it has prob > p_thresh and overlaps the node, then the node has support
	    if (kmer_path_probs[i]>p_thresh)
	    {
	        for (deque<Interval>::const_iterator it=prgs[pnode->id]->kmer_paths[i].path.begin(); it!=prgs[pnode->id]->kmer_paths[i].path.end(); ++it)
                {
		    if (it->end > n->second->pos.start and it->start < n->second->pos.end)
		    {
			prgs[pnode->id]->prg.add_read_support_node(n->second);
			found_support = true;
			num_supported_nodes += 1;
			break; //breaks loop over intervals in minimizer
		    }
                }
	    }
	    ++i;
	}
    }
    cout << "found " << num_supported_nodes << " supported nodes" << endl;
    assert(num_supported_nodes>0);
    assert(num_supported_nodes<kmer_path_hit_counts.size());

    // now infer the paths through these supported nodes
    // for now, write them straight to a fasta
    // in future, could choose most likely of these paths
    cout << "next infer the read supported graph (connecting the supported nodes" << endl;
    prgs[pnode->id]->prg.infer_read_supported_graph();
    cout << "finally write the paths through the supported graph to fasta" << endl;
    prgs[pnode->id]->prg.write_read_supported_graph_paths("supported_path_seqs" + prgs[pnode->id]->name + ".fasta");
    return;
}

uint32_t nchoosek (uint32_t n, uint32_t k)
{
    assert(n >= k);
    //assert(n >= 0);
    //assert(k >= 0);

    if (n == 0) {return 0;}
    if (k == 0) {return 1;}
    if (n == k) {return 1;}
    
    return (n * nchoosek(n - 1, k - 1)) / k;
}
