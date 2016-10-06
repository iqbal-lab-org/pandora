#ifndef __SEQ_H_INCLUDED__   // if seq.h hasn't been included yet...
#define __SEQ_H_INCLUDED__

#include <string>
#include <vector>
#include <ostream>
#include "minimizer.h"

using std::vector;
using namespace std;

class Seq {
  public:
    uint32_t id;
    string name;
    string seq;
    //vector<Minimizer*> sketch; //argument for set - removes duplicates and orders them
    set<Minimizer*, pMiniComp> sketch;
    Seq(uint32_t, string, string, uint32_t, uint32_t);
    ~Seq();
    void minimizer_sketch (uint32_t w, uint32_t k);
  friend ostream& operator<< (ostream& out, const Seq& data); 
};

#endif
