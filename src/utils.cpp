#include <cstring>
#include <vector>
#include <iostream>
//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include "utils.h"
#include "index.h"
#include "localPRG.h"
#include "seq.h"
#include "pangraph.h"
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

void infer_localPRG_order_for_reads(MinimizerHits* minimizer_hits, PanGraph* pangraph, const int max_diff, const uint32_t cluster_thresh, const uint32_t k)
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
        //cout << "Hit: " << **mh_previous << endl;
        if((*mh_current)->read_id!=(*mh_previous)->read_id or (*mh_current)->prg_id!=(*mh_previous)->prg_id or (*mh_current)->strand!=(*mh_previous)->strand or (abs((int)(*mh_current)->read_interval.start - (int)(*mh_previous)->read_interval.start)) > max_diff)
        {
            if (current_cluster.size() > cluster_thresh)
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
    if (current_cluster.size() > cluster_thresh)
    {
        clusters_of_hits.insert(current_cluster);
        //cout << "Found final cluster of size: " << current_cluster.size() << " ending with " << **mh_previous << endl;
    } else {
        //cout << "Final cluster not added had size " << current_cluster.size() << " and ended with " << **mh_previous << endl;
    }

    // Next order clusters, remove contained ones, and add inferred order to pangraph
    if (clusters_of_hits.size() == 0) { return;}
    set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_previous = clusters_of_hits.begin();
    pangraph->add_node((*(*c_previous).begin())->prg_id, (*(*c_previous).begin())->read_id);
    for (set<set<MinimizerHit*, pComp>, clusterComp>::iterator c_current = ++clusters_of_hits.begin(); c_current != clusters_of_hits.end(); ++c_current)
    {
        if(((*(*c_current).begin())->read_id == (*(*c_previous).begin())->read_id) &&  ((*(*c_current).begin())->prg_id != (*(*c_previous).begin())->prg_id) && ((*--(*c_current).end())->read_interval.start > (*--(*c_previous).end())->read_interval.start + k - 1) ) // NB we expect noise in the k-1 kmers overlapping the boundary of two clusters, so force the next cluster to have at least a hit which is outside this region
        {
            pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id);
	    pangraph->add_edge((*(*c_previous).begin())->prg_id, (*(*c_current).begin())->prg_id);
            c_previous = c_current;
        //} else {
        //    cout << "Contained cluster not added to order" << endl;
        } else if ((*(*c_current).begin())->read_id != (*(*c_previous).begin())->read_id)
	{
	    // if we just started looking at hits for a new read, add the first cluster
	    pangraph->add_node((*(*c_current).begin())->prg_id, (*(*c_current).begin())->read_id);
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
		    infer_localPRG_order_for_reads(mh, pangraph, max_diff, cluster_thresh, k);
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
            infer_localPRG_order_for_read(mh, pangraph, max_diff, cluster_thresh, k);
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
        infer_localPRG_order_for_reads(mh, pangraph, max_diff, cluster_thresh, k);
        now = time(0);
        dt = ctime(&now);
        sdt = dt.substr(0,dt.length()-1);
        cout << sdt << " Update coverages within genes" << endl;
        update_covgs_from_hits(prgs, mh);
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


/*void infer_most_likely_prg_path_for_read(const vector<LocalPRG*>& prgs, MinimizerHits* mhs, 
*/
