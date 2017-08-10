#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "pannode.h"
//#include "seq.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const string n): id(i), name(n) 
{
    vector<PanNode*> v;
    unordered_map<uint32_t,uint16_t> m;
    for (uint k=0; k<4; ++k)
    {
        outNodes.push_back(v);
	outNodeCounts.push_back(m);
    }
}

void PanNode::add_read(const uint32_t j)
{
    foundReads.push_back(j);
}

void PanNode::add_hits(const set<MinimizerHit*, pComp>& c)
{
    //cout << "adding " << c.size() << " hits" << endl;
    foundHits.insert(c.begin(), c.end());
    //cout << "added " << foundHits.size() << " hits" << endl;
}

bool PanNode::operator == (const PanNode& y) const {
    if (id!= y.id) {return false;}

    bool found_outnode;
    for (uint k=0; k!=3; ++k)
    {
    	for (uint32_t i=0; i!=outNodes[k].size(); ++i)
    	{
	    found_outnode = false;	
	    for (uint32_t j=0; j!=y.outNodes[k].size(); ++j)
	    {
	    	if (outNodes[k][i]->id == y.outNodes[k][j]->id)
	    	{
		    found_outnode = true;
		    break;
	    	}
	    }
	    if (found_outnode == false)
	    {
	        return false;
	    }	
        }
        for (uint32_t i=0; i!=y.outNodes[k].size(); ++i)
        {   
            found_outnode = false;  
            for (uint32_t j=0; j!=outNodes[k].size(); ++j)
            {   
                if (y.outNodes[k][i]->id == outNodes[k][j]->id)
                {   
                    found_outnode = true;
                    break;
                }
            }
            if (found_outnode == false)
            {   
                return false;
            }   
        }
    }
    return true;
}

std::ostream& operator<< (std::ostream & out, PanNode const& m) {
    out << m.id << " -> ";
    for (uint32_t i=0; i!=m.outNodes[3].size(); ++i)
    { out << m.outNodes[3][i]->id << "++ ";}
    for (uint32_t i=0; i!=m.outNodes[1].size(); ++i)
    { out << m.outNodes[1][i]->id << "+- ";}
    for (uint32_t i=0; i!=m.outNodes[2].size(); ++i)
    { out << m.outNodes[2][i]->id << "-+ ";}
    for (uint32_t i=0; i!=m.outNodes[0].size(); ++i)
    { out << m.outNodes[0][i]->id << "-- ";}
    return out ;
}
