#ifndef __SEQ_H_INCLUDED__   // if seq.h hasn't been included yet...
#define __SEQ_H_INCLUDED__

#include <string>
#include <set>
//#include<functional>
#include "minimizer.h"

using std::set;
using namespace std;

class Seq {
  public:
    uint32_t id;
    string name;
    string seq;
    set<Minimizer*, pMiniComp> sketch; //argument for set - removes duplicates and orders them
    Seq(uint32_t, string, string, uint32_t, uint32_t);
    ~Seq();
    void minimizer_sketch (uint32_t w, uint32_t k);
    void print() const;
};

#endif
