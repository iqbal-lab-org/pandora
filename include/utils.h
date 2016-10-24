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

void index_prg_file(vector<LocalPRG*>&, string, Index*, uint32_t, uint32_t);
void add_read_hits(uint32_t, string, string, MinimizerHits*, Index*, uint32_t, uint32_t);
void infer_localPRG_order_for_read(MinimizerHits*, PanGraph*, int, uint32_t, uint32_t);
void pangraph_from_read_file(string, PanGraph*, Index*, uint32_t, uint32_t, int, uint32_t);

#endif
