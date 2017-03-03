#include <iostream>
#include <string>
#include <algorithm>
#include "kmernode.h"
#include "path.h"
#include "utils.h" // for pointer_values_equal

using namespace std;

KmerNode::KmerNode (uint32_t i, const Path& p): id(i), path(p), covg({0,0}) {}

std::ostream& operator<< (std::ostream & out, KmerNode const& n) {
    out << n.id << " " << n.path << endl;
    for (uint32_t i=0; i!=n.outNodes.size(); ++i)
    {
        out << n.id << " -> " << n.outNodes[i]->id << endl;
    }
    return out ;
}

bool KmerNode::operator == (const KmerNode& y) const {
    if (!(path == y.path)) {return false;}
/*    if (outNodes.size() != y.outNodes.size()) {return false;}
    if (inNodes.size() != y.inNodes.size()) {return false;}
    for (uint32_t i=0; i!=outNodes.size(); ++i)
    {
        pointer_values_equal<KmerNode> eq = { outNodes[i] };
        if ( find_if(y.outNodes.begin(), y.outNodes.end(), eq) == y.outNodes.end() )
	{return false;}
    }
    for (uint32_t i=0; i!=inNodes.size(); ++i)
    {
        pointer_values_equal<KmerNode> eq = { inNodes[i] };
        if ( find_if(y.inNodes.begin(), y.inNodes.end(), eq) == y.inNodes.end() )
        {return false;}
    }*/
    return true;
}
