/*
 * Adapted from
 * https://github.com/iqbal-lab/cortex/blob/master/include/cortex_var/many_colours/element.h
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
  element.h  - defines the interface for the de Bruijn graph node. The
  implementation is complemented by a hash table that stores every node indexed
  by kmers (BinaryKmers).

  The element routines, ie the one required by hash_table/priority queue, are
  prefixed with element_

  The de Bruijn based routines are prefixed with db_node_
*/
#ifndef CORTEX_DB_NODE_H
#define CORTEX_DB_NODE_H

#include "ns.h"
#include "cortex/binary_kmer.h"

using Covg = uint32_t;
using Edges = char;

class cortex::Node {
    // todo: figure out private/public
    const BinaryKmer kmer;
    const Covg coverage;
    const Edges individual_edges;
    const char status; // will cast a NodeStatus to char
    const char allele_status;
};

#endif // CORTEX_DB_NODE_H
