#include <iostream>
#include <string>
#include "msarecord.h"

using namespace std;

MSARecord::MSARecord (uint32_t i, string n, string p): id(i), name(n), seq(p) {
}

std::ostream& operator<< (std::ostream & out, MSARecord const& data) {
    out << data.name;
    return out ;
}

