#ifndef CORTEX_NUCLEOTIDE_H
#define CORTEX_NUCLEOTIDE_H

#include "cortex/ns.h"

using cortex::Nucleotide;

enum class cortex::Nucleotide {
    Adenine,
    Cytosine,
    Guanine,
    Thymine,
    Undefined,
};

// char to binary nucleotide
Nucleotide cton(char c);

// binary nucleotide to char
char ntoc(Nucleotide n);

// nucleotide complement
Nucleotide complement(Nucleotide n);

// boolean good_symbol(char c);

#endif // CORTEX_NUCLEOTIDE_H
