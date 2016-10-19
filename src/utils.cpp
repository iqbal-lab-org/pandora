#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>

#include "utils.h"
#include "index.h"
#include "localPRG.h"
#include "seq.h"
#include "pangraph.h"
#include "minihits.h"
#include "minihit.h"

void index_prg_file(vector<LocalPRG*>& prgs, string filepath, Index* idx, uint32_t w, uint32_t k)
{
    uint32_t id = 0;
    string name, read, line;
    LocalPRG *s;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!name.empty() && !read.empty())
                {
                    //cout << name << ": " << read << endl;
                    s = new LocalPRG(id, name, read);
		    s->minimizer_sketch(idx, w, k);
                    
                    if (s!=nullptr)
                        {prgs.push_back(s);}
                }
                name.clear();
                read.clear();
		id++;
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
        if (!name.empty() && !read.empty())
        {
            //cout << name << ": " << read << endl;
            s = new LocalPRG(id, name, read);
            s->minimizer_sketch(idx, w, k);
            if (s!=nullptr)
                {prgs.push_back(s);}
        }
        cout << "Number of LocalPRGs added to Index: " << prgs.size() << endl;
        myfile.close();
    } else {
        cout << "Unable to open PRG file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

void add_read_hits(uint32_t id, string name, string seq, MinimizerHits* hits, Index* idx, uint32_t w, uint32_t k)
{
    // creates Seq object for the read, then looks up minimizers in the Seq sketch and adds hits to a global MinimizerHits object
    Seq s = Seq(id, name, seq, w, k);
    for(set<Minimizer*, pMiniComp>::iterator it = s.sketch.begin(); it != s.sketch.end(); ++it)
    {
        if (idx->minhash.find((*it)->kmer) != idx->minhash.end())
        {
	    for (vector<MiniRecord>::iterator it2=idx->minhash[(*it)->kmer].begin(); it2!=idx->minhash[(*it)->kmer].end(); ++it2)
            {
	        hits->add_hit(s.id, *it, *it2, 1);
            }
        }
    }
    return;
}

void infer_localPRG_order_for_read(MinimizerHits* minimizer_hits, PanGraph* pangraph, int max_diff, uint32_t cluster_thresh)
{
    // this step infers the gene order for a read
    // orders hits from a set of minimizer hits, clusters them, removes noise hits, and adds the inferred gene order to the pangraph
    set<set<MinimizerHit*, pComp>,clusterComp> clusters_of_hits;

    if (minimizer_hits->hits.size() == 0) {return;}

    // First cluster hits matching same localPRG, not more than max_diff read bases from the last hit. 
    // Keep clusters with more than cluster_thresh hits.
    set<MinimizerHit*, pComp>::iterator mh_previous = minimizer_hits->hits.begin();
    set<MinimizerHit*, pComp> current_cluster;
    current_cluster.insert(*mh_previous);
    for (set<MinimizerHit*, pComp>::iterator mh_current = ++minimizer_hits->hits.begin(); mh_current != minimizer_hits->hits.end(); ++mh_current)
    {
        if((*mh_current)->prg_id!=(*mh_previous)->prg_id or ((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start) > max_diff)
        {
            if (current_cluster.size() > cluster_thresh)
            {
                //cout << "Found cluster of size: " << current_cluster.size() << endl;
                clusters_of_hits.insert(current_cluster);
            }
            current_cluster.clear();
            current_cluster.insert(*mh_current);
        } else {
            current_cluster.insert(*mh_current);
        }
        mh_previous = mh_current;
    }
    if (current_cluster.size() > cluster_thresh)
    {
        clusters_of_hits.insert(current_cluster);
        //cout << "Found cluster of size: " << current_cluster.size() << endl;
    }

    // Next order clusters, remove contained ones, and add inferred order to pangraph
    if (clusters_of_hits.size() == 0) { return;}
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id);
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start)
        {
            pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_previous).begin())->read_id);
            pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        //} else {
            //cout << "Contained cluster not added to order" << endl;
        }
    }
    return;
}

void pangraph_from_read_file(string filepath, PanGraph* pangraph, Index* idx, uint32_t w, uint32_t k)
{
    string name, read, line;
    uint32_t id = 0;
    MinimizerHits* mh;

    ifstream myfile (filepath);
    if (myfile.is_open())
    {
        while ( getline (myfile,line).good() )
        {
            if (line.empty() || line[0] == '>' )
            {
                if (!name.empty() && !read.empty())
                {
                    mh = new MinimizerHits;
                    add_read_hits(id, name, read, mh, idx, w, k);
		    infer_localPRG_order_for_read(mh, pangraph, 500, 4);
                }
                name.clear();
                read.clear();
		id++;
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
        if (!name.empty() && !read.empty())
        {
            mh = new MinimizerHits;
            add_read_hits(id, name, read, mh, idx, w, k);
            infer_localPRG_order_for_read(mh, pangraph, 500, 4);
        }
        //cout << "Number of reads found: " << id+1 << endl;
        myfile.close();
    } else {
        cerr << "Unable to open read file " << filepath << endl;
        exit (EXIT_FAILURE);
    }
    return;
}

