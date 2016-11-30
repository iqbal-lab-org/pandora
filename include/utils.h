#ifndef __UTILS_H_INCLUDED__   // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

#include <vector>
#include "localPRG.h"

class Index;
class MinimizerHits;
class PanGraph;

template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};

void index_prg_file(vector<LocalPRG*>&, const string&, Index*, const uint32_t, const uint32_t);
void add_read_hits(uint32_t, const string&, const string&, MinimizerHits*, Index*, const uint32_t, const uint32_t);
void infer_localPRG_order_for_reads(MinimizerHits*, PanGraph*, const int, const uint32_t, const uint32_t);
void pangraph_from_read_file(const string&, PanGraph*, Index*, const vector<LocalPRG*>&, const uint32_t, const uint32_t, const int, const uint32_t);
void update_covgs_from_hits(const vector<LocalPRG*>&, MinimizerHits*);

#endif
