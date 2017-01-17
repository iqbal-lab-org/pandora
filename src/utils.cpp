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

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

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
    cout << sdt << " Found " << hit_count << " hits found for read " << name << " so size of MinimizerHits is now " << hits->hits.size() << endl;
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
	    //cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
            if (pn < 0.001)
            {
                cout << "Found cluster of size: " << current_cluster.size() << " for prg " << (*mh_previous)->prg_id << endl;
                clusters_of_hits.insert(current_cluster);
		/*for (set<MinimizerHit*, pComp>::iterator mh = current_cluster.begin(); mh!=current_cluster.end(); ++mh)
		{
		    cout << **mh << endl;
		}
		cout << endl;*/
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
    //cout << "pnull is " << pn << " for cluster of size " << current_cluster.size() << endl;
    if (pn < 0.001)
    {
        clusters_of_hits.insert(current_cluster);
        cout << "Found final cluster of size: " << current_cluster.size() << " for prg " << (*mh_previous)->prg_id << endl;
    } else {
        //cout << "Final cluster not added had size " << current_cluster.size() << " and ended with " << **mh_previous << endl;
    }
    cout << "Found " << clusters_of_hits.size() << " clusters" << endl;

    // Next order clusters, remove contained ones, and add inferred order to pangraph    
    if (clusters_of_hits.size() == 0) { return;}
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id, *c_previous);
    cout << "first cluster added " << (*(*c_previous).begin())->prg_id << endl; 
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) &&  ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) && ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start) ) // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, so force the next cluster to have at least a hit which is outside this region
        {
	    cout << "added cluster with id " << (*(*c_current).begin())->prg_id << endl;
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
	} else {
	    cout << "Did not add cluster. Criteria which may have failed:" << endl;
	    cout << "read_ids equal: " << (*(*c_current).begin())->read_id << "==" << (*(*c_previous).begin())->read_id << ", " << ((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) << endl;
	    cout << "prg_ids different: " << (*(*c_current).begin())->prg_id << "!=" << (*(*c_previous).begin())->prg_id << ", " << ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) << endl;
	    cout << "not contained: " << (*--(*c_current).end())->read_interval.start << ">" << (*--(*c_previous).end())->read_interval.start << ", " << ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start) << endl;
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
                    //mh = new MinimizerHits();
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
            //mh = new MinimizerHits();
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
	    prgs[pnode->second->id]->infer_most_likely_prg_paths_for_corresponding_pannode(pnode->second, k, 0.00001);
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
        for (uint32_t i=0; i!= prgs.size(); ++i)
	{
	    if (prgs[i]->id == (*mh)->prg_id)
	    {
	        prgs[i]->update_covg_with_hit(*mh);
		break;
	    }
	}
	
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
            cout << "p_null = " << p << ", using |y|=" << prgs[i]->kmer_paths.size() << " and |x|=" << cluster_of_hits.size() << " for prg " << prgs[i]->id << " and strand " << (*cluster_of_hits.begin())->strand << endl;
            return p;
	}
    }
    // if we have got here, then it couldn't be found
    cout << "did not find id" << endl;
    assert(id_found == true);
    return p;
}

/*void infer_most_likely_prg_path_for_pannode(const vector<LocalPRG*>& prgs, PanNode* pnode, uint32_t k, float e_rate)
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

    //now we iterate through the graph from the outmost level to the lowest level working out the most likely path(s)
    //need 2 data structures, one to remember what the most probable path(s) were for var sites already considered
    map<uint32_t, vector<pair<vector<LocalNode*>, float>>> max_path_index; // maps from pre_site_ids seen before to a vector of pairs representing max node_paths and their probs (should be the same)
    //and the second remembering which minimizing kmers have been seen/used before so not included multiple times in the probability
    vector<bool> kmer_considered_before(prgs[pnode->id]->kmer_paths.size(),false);

    // start with the outmost level
    uint8_t max_level = 0;
    for (auto const& element : prgs[pnode->id]->prg.index) {
        max_level = max(max_level, element.first);
    }
    // and for each level..
    for (uint8_t level = max_level; level < max_level + 1; --level)
    {
        // ...for each varsite at this level...
        for (uint i = 0; i!=prgs[pnode->id]->prg.index[level].size(); ++i)
        {
	    // ...find the maximally probable paths through varsite
            vector<pair<vector<LocalNode*>, float>> u, v, w;
	    vector<LocalNode*> x;

            // add the first node of each alternate allele for the varsite to a vector
            uint32_t pre_site_id = prgs[pnode->id]->prg.index[level][i].first, post_site_id = prgs[pnode->id]->prg.index[level][i].second;
            for (uint j = 0; j!=prgs[pnode->id]->prg.nodes[pre_site_id]->outNodes.size(); ++j)
	    {
		x.push_back(prgs[pnode->id]->prg.nodes[pre_site_id]->outNodes[j]);
		u.push_back(make_pair(x, 1));
		x.clear();
	    } 

            // then until we reach the end varsite:
	    while (u.size()>0)
	    {
                // for each pair in u
                for (uint j = 0; j!=u.size(); ++j)
		{
            	    // if the path in the pair ends at the end of the varsite, it is a done path, add to w
            	    if (u[j].first.back()->id == post_site_id)
		    {
			w.push_back(u[j]);
		    } else {
            		// otherwise look up the last node in the max_path_index
            		// if it is there, then we can multiple the running total prob for this allele, 
            		// and extend the path for each of the maximal paths through the sub_varsite
            		// and add the updated path/prob pairs to v
            		map<uint32_t, vector<pair<vector<LocalNode*>, float>>>::iterator it=max_path_index.find(u[j].first.back()->id);
            		if (it != max_path_index.end())
			{
			    //extend node path with the max paths seen before from this point
			    for (uint n = 0; n!=it->second.size(); ++n)
			    {
				v.push_back(u[j]);
				v.back().first.insert(v.back().first.end(), it->second[n].first.begin(), it->second[n].first.end());
				v.back().second = v.back().second * it->second[n].second;
			    }
			} else {
            		    // if not, then work out which minhits overlap this node, and multiply their probabilities, 
            		    // then add the outnode to the end (should only be 1) to get a new path/prob pair to add to v
                            // during this process, update kmer_considered_before to reflect the minihits now used in probabilities
                        }
		    }
		}
		// once done for all of what was in u, set u = v
		u = v;
	    }
	    // when u empty, should have final set in w and can work out the max of the probs
	    float max_prob = 0;
	    // find max for all probs we could work out values for
    	    for (uint n = 0; n!=w.size(); ++n)
	    {
		if (w[n].second != 1) // if no mini kmers in prg overlapping node, can't define prob data came from that node, 
			    		           // so set the prob to 1 intially, then set to max path value 
						   // only happens for really short paths
		{
		    max_prob = max(max_prob, w[n].second);
		}
	    }
            // now add both those paths achieving max, and paths cannot judge to max_path_index
	    max_path_index[pre_site_id] = u; // know u is empty vector
            for (uint n = 0; n!=w.size(); ++n)
            {
		if (w[n].second == 1 or w[n].second == max_prob)
		{
		    max_path_index[pre_site_id].push_back(w[n]); // know u is empty vector
		}
	    }
	}
    }
    
    // finally, write out the maximal paths to file
    write_paths_to_fasta(max_path_index[0], "supported_path_seqs_" + prgs[pnode->id]->name + ".fasta");

    return;
}*/

uint32_t nchoosek (uint32_t n, uint32_t k)
{
    assert(n >= k || assert_msg("Currently the model assumes that the most a given kmer (defined by position) can occur is once per read, i.e. an error somewhere else in the read cannot result in this kmer. If you are getting this message, then you have evidence of violation of this assumption. Either try using a bigger k, or come up with a better model"));
    //assert(n >= 0);
    //assert(k >= 0);

    if (n == 0) {return 0;}
    if (k == 0) {return 1;}
    if (n == k) {return 1;}
    
    return (n * nchoosek(n - 1, k - 1)) / k;
}

/*void write_paths_to_fasta(vector<pair<vector<LocalNode*>, float>> p, const string& filepath)
{
    ofstream handle;
    handle.open (filepath);
    for (uint i = 0; i!= p.size(); ++i)
    {
        handle << ">\t" << it->second->id << "\t";
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
}*/
