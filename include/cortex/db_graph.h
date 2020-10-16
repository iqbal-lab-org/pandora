/*
 * Adapted from
 * https://github.com/iqbal-lab/cortex/blob/master/include/cortex_var/core/dB_graph.h
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
  dB_graph.h

  all the routines as prefixed with db_graph
*/

#ifndef DB_GRAPH_H_
#define DB_GRAPH_H_

#include <stdio.h>
#include <unordered_map>
#include "cortex/ns.h"
#include "cortex/binary_kmer.h"
#include "cortex/db_node.h"

class cortex::DBGraph {
private:
    std::unordered_map<cortex::BinaryKmer, cortex::Node> hash_table;

public:
    //    int db_graph_clip_tip(cortex::Node* node, int limit, cortex::DBGraph*
    //    db_graph);
    int clip_tips(int limit);

    //    void db_graph_set_all_visited_nodes_to_status_none(cortex::DBGraph*
    //    hash_table);
    void set_all_visited_to_none();

    // pays no attention to whether there is an edge joining current_node to the node
    // you would get by adding this nucleotide. just checks to see if such a node is in
    // the graph
    //    dBNode* db_graph_get_next_node(dBNode* current_node,
    //        Orientation current_orientation, Orientation* next_orientation, Nucleotide
    //        edge, Nucleotide* reverse_edge, dBGraph* db_graph);
    cortex::Node& get_next_node(/*todo*/);
};

#endif /* DB_GRAPH_H_ */
