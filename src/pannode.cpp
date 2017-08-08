#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "pannode.h"
//#include "seq.h"
#include "utils.h"

using namespace std;

PanNode::PanNode (const uint32_t i, const string n): id(i), name(n) {}

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
    for (uint32_t i=0; i!=outNodes.size(); ++i)
    {
	found_outnode = false;	
	for (uint32_t j=0; j!=y.outNodes.size(); ++j)
	{
	    if (outNodes[i]->id == y.outNodes[j]->id)
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
    for (uint32_t i=0; i!=y.outNodes.size(); ++i)
    {   
        found_outnode = false;  
        for (uint32_t j=0; j!=outNodes.size(); ++j)
        {   
            if (y.outNodes[i]->id == outNodes[j]->id)
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
    return true;
}

std::ostream& operator<< (std::ostream & out, PanNode const& m) {
    out << m.id << " -> ";
    for (uint32_t i=0; i!=m.outNodes.size(); ++i)
    { out << m.outNodes[i]->id << " ";}
    return out ;
}
