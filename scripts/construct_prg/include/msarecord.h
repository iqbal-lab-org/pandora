#ifndef __MSARECORD_H_INCLUDED__   // if msarecord.h hasn't been included yet...
#define __MSARECORD_H_INCLUDED__

#include <string>
#include <ostream>
#include "msarecord.h"

using namespace std;

class MSARecord {
  public:
    uint32_t id;
    string name;
    string seq;
  friend ostream& operator<< (ostream& out, const MSARecord& data); 
};

#endif
