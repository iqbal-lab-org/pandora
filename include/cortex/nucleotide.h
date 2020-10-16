/*
 * Adapted from https://github.com/iqbal-lab/cortex/blob/master/include/basic/base_encoding/event_encoding.h
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 *
 * CORTEX project contacts:
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  base_encoding/event_encoding.h
*/
#ifndef CORTEX_NUCLEOTIDE_H
#define CORTEX_NUCLEOTIDE_H

#include <string>
#include <vector>
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

// nucleotides to string
std::string ntos(std::vector<Nucleotide>& nucleotides);

// nucleotide complement
Nucleotide complement(Nucleotide n);

// boolean good_symbol(char c);

#endif // CORTEX_NUCLEOTIDE_H
