#include "localgraph.h"

LocalGraph::LocalGraph() { }

LocalGraph::~LocalGraph() { nodes.clear(); }
// add the node with a given id and seq if the id is not already in the nodes
void LocalGraph::add_node(
    const uint32_t& id, const std::string& seq, const Interval& pos)
{
    const bool sequence_and_interval_length_match = seq.length() == pos.length;
    if (!sequence_and_interval_length_match) {
        fatal_error("Error adding node to Local Graph: sequence and interval length do "
                    "not match");
    }

    auto it = nodes.find(id);
    if (it == nodes.end()) {
        LocalNodePtr n(std::make_shared<LocalNode>(seq, pos, id));
        nodes[id] = n; // add the node to the map
        // add the node to the interval indexes for fast overlap queries
        if (pos.length == 0)
            startIndexOfZeroLengthIntervals[pos.start] = n;
        else
            intervalTree.add(pos.start, pos.get_end(), n);
        startIndexOfAllIntervals[pos.start] = n;
    } else {
        const bool node_with_same_id_seq_and_pos_already_added
            = (it->second->seq == seq) && (it->second->pos == pos);
        if (!node_with_same_id_seq_and_pos_already_added) {
            fatal_error("Error adding node to Local Graph: node with ID ", id,
                " already exists in graph, but with different sequence or pos");
        }
    }
}

void LocalGraph::add_edge(const uint32_t& from, const uint32_t& to)
{
    auto from_it = nodes.find(from);
    auto to_it = nodes.find(to);

    const bool both_nodes_exist = (from_it != nodes.end()) && (to_it != nodes.end());
    if (!both_nodes_exist) {
        fatal_error("Cannot add edge to Local Graph: source (", from, ") or target (",
            to, ") node does not exist");
    }

    LocalNodePtr f = (nodes.find(from)->second);
    LocalNodePtr t = (nodes.find(to)->second);

    const bool nodes_do_not_overlap = f->pos.get_end() > t->pos.start;
    if (nodes_do_not_overlap) {
        fatal_error(
            "Cannot add edge to Local Graph: source and target nodes do not overlap");
    }
    f->outNodes.push_back(t);
}

void LocalGraph::write_gfa(const std::string& filepath) const
{
    std::ofstream handle;
    handle.open(filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << std::endl;
    for (const auto& node : nodes) {
        handle << "S\t" << node.second->id << "\t";
        if (node.second->seq.empty()) {
            handle << "*";
        } else {
            handle << node.second->seq;
        }
        handle << "\tRC:i:" << node.second->covg << std::endl;
        for (uint32_t j = 0; j < node.second->outNodes.size(); ++j) {
            handle << "L\t" << node.second->id << "\t+\t"
                   << node.second->outNodes[j]->id << "\t+\t0M" << std::endl;
        }
    }
    handle.close();
}

void LocalGraph::read_gfa(const std::string& filepath)
{
    uint32_t id, from, to;
    std::string line;
    std::vector<std::string> split_line;
    uint32_t i = 0;

    std::ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");

                const bool line_is_consistent = split_line.size() >= 3;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                if (split_line[2] == "*") {
                    split_line[2] = "";
                }
                id = std::stoi(split_line[1]);
                add_node(id, (std::string)split_line[2],
                    Interval(i, i + split_line[2].size()));
                i += split_line[2].size();
            }
        }

        myfile.clear();
        myfile.seekg(0, myfile.beg);
        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");

                const bool line_is_consistent = split_line.size() >= 5;
                if (!line_is_consistent) {
                    fatal_error("Error reading GFA. Offending line: ", line);
                }

                if (split_line[2] == split_line[4]) {
                    from = stoi(split_line[1]);
                    to = stoi(split_line[3]);
                } else {
                    from = stoi(split_line[3]);
                    to = stoi(split_line[1]);
                }
                add_edge(from, to);
            }
        }
    } else {
        fatal_error("Unable to open GFA file: ", filepath);
    }
}

std::vector<PathPtr> LocalGraph::walk(
    const uint32_t& node_id, const uint32_t& pos, const uint32_t& len) const
{ // node_id: where to start the walk, pos: the position in the node_id, len = k+w-1 ->
  // the length that the walk has to go through - we are sketching kmers in a graph
    // walks from position pos in node node for length len bases
    const bool pos_exists_in_node = (nodes.at(node_id)->pos.start <= pos)
        && (pos <= nodes.at(node_id)->pos.get_end());
    if (!pos_exists_in_node) {
        fatal_error("Error walking Local Graph: pos ", pos, " does not exist in node ",
            node_id);
    }

    std::vector<PathPtr> return_paths, walk_paths;
    return_paths.reserve(20);
    walk_paths.reserve(20);
    prg::Path p, p2;
    std::deque<Interval> d;

    if (pos + len
        <= nodes.at(node_id)
               ->pos.get_end()) { // checks if we can go until the end of the kmer
        p.initialize(Interval(
            pos, pos + len)); // create the path containing this interval of the node
        return_paths.push_back(std::make_shared<prg::Path>(p));
        return return_paths;
    }
    uint32_t len_added = std::min(nodes.at(node_id)->pos.get_end() - pos, len);

    if (len_added < len) {
        for (auto it = nodes.at(node_id)->outNodes.begin();
             it != nodes.at(node_id)->outNodes.end(); ++it) {
            walk_paths = walk((*it)->id, (*it)->pos.start, len - len_added);
            for (auto& walk_path : walk_paths) {
                // Note, would have just added start interval to each item in
                // walk_paths, but can't seem to force result of it2 to be non-const
                p2.initialize(Interval(pos, nodes.at(node_id)->pos.get_end()));
                p2.insert_to_the_end(walk_path->begin(), walk_path->end());
                if (p2.length() == len) {
                    return_paths.push_back(std::make_shared<prg::Path>(p2));
                }
            }
        }
    }
    return return_paths;
}

std::vector<PathPtr> LocalGraph::walk_back(
    const uint32_t& node_id, const uint32_t& pos, const uint32_t& len) const
{
    // walks from position pos in node back through prg for length len bases
    const bool pos_exists_in_node = (nodes.at(node_id)->pos.start <= pos)
        && (pos <= nodes.at(node_id)->pos.get_end());
    if (!pos_exists_in_node) {
        fatal_error("Error walking Local Graph: pos ", pos, " does not exist in node ",
            node_id);
    }

    std::vector<PathPtr> return_paths, walk_paths;
    return_paths.reserve(20);
    walk_paths.reserve(20);
    prg::Path p, p2;
    std::deque<Interval> d;

    if (nodes.at(node_id)->pos.start + len <= pos) {
        p.initialize(Interval(pos - len, pos));
        return_paths.push_back(std::make_shared<prg::Path>(p));
        return return_paths;
    }

    uint32_t len_added = std::min(pos - nodes.at(node_id)->pos.start, len);

    std::vector<LocalNodePtr>::iterator innode;
    if (len_added < len) {
        for (auto it = nodes.begin(); it != nodes.find(node_id); ++it) {
            innode = find(it->second->outNodes.begin(), it->second->outNodes.end(),
                nodes.at(node_id));
            if (innode != it->second->outNodes.end()) {
                walk_paths = walk_back(
                    it->second->id, it->second->pos.get_end(), len - len_added);
                for (uint32_t i = 0; i != walk_paths.size(); ++i) {
                    p2.initialize(*(walk_paths[i]));
                    p2.add_end_interval(Interval(nodes.at(node_id)->pos.start, pos));
                    if (p2.length() == len) {
                        return_paths.push_back(std::make_shared<prg::Path>(p2));
                    }
                }
            }
        }
    }
    return return_paths;
}

LocalNodePtr LocalGraph::get_previous_node(const LocalNodePtr n) const
{
    // returns a previous node if there is one/many
    if (n->id == 0) {
        return nullptr;
    } else {
        for (const auto& c : nodes) {
            if (find(c.second->outNodes.begin(), c.second->outNodes.end(), n)
                != c.second->outNodes.end()) {
                return c.second;
            } else if (c.first > n->id) {
                break;
            }
        }
        // if we get here, there was no previous node to be found.
        return nullptr;
    }
}

std::vector<LocalNodePtr> LocalGraph::nodes_along_string(
    const std::string& query_string, bool end_to_end) const
{
    // Note expects the query string to start at the start of the PRG - can change this
    // later
    const bool graph_is_empty = nodes.empty();
    if (graph_is_empty) {
        fatal_error("Error getting nodes along a sequence: graph is empty");
    }

    std::vector<std::vector<LocalNodePtr>> u, v, w; // u <=> v -> w
    // ie reject paths in u, or extend and add to v
    // then set u=v and continue
    // final output w, which is filtered and a path returned
    u.reserve(100);
    v.reserve(100);
    w.reserve(100);
    std::vector<LocalNodePtr> npath;
    std::string candidate_string = "";
    bool extended = true;

    // if there is only one node in PRG, simple case, do simple string compare
    if (nodes.size() == 1
        and strcasecmp(query_string.c_str(), nodes.at(0)->seq.c_str()) == 0) {
        return { nodes.at(0) };
    }

    u = { { nodes.at(0) } };

    while (!u.empty()) {
        for (const auto& p : u) {
            candidate_string = "";
            for (const auto& s : p) {
                candidate_string += s->seq;
            }

            for (uint32_t j = 0; j != p.back()->outNodes.size(); ++j) {
                // if the start of query_string matches extended candidate_string, want
                // to query candidate path extensions
                auto comp_string = candidate_string + p.back()->outNodes[j]->seq;
                auto comp_length = std::min(query_string.size(), comp_string.size());
                if (strcasecmp(query_string.substr(0, comp_length).c_str(),
                        comp_string.substr(0, comp_length).c_str())
                    == 0) {
                    if ((!end_to_end
                            and candidate_string.size()
                                    + p.back()->outNodes[j]->seq.size()
                                >= query_string.size())
                        or p.back()->outNodes[j]->outNodes.empty()) {
                        // we have now found the whole of the query_string or reached
                        // end of graph
                        auto p_copy = p;
                        p_copy.push_back(p_copy.back()->outNodes[j]);
                        while (!p_copy.back()->outNodes.empty() and extended) {
                            extended = false;
                            for (const auto& n : p_copy.back()->outNodes) {
                                if (n->pos.length == 0) {
                                    p_copy.push_back(n);
                                    extended = true;
                                    break;
                                }
                            }
                        }
                        w.push_back(p_copy);
                        continue;
                    } else {
                        v.push_back(p);
                        v.back().push_back(p.back()->outNodes[j]);
                    }
                }
            }
        }
        u = v;
        // sanity check, we don't care enough about finding the path to have exploding
        // vectors
        if (u.size() > 10000)
            u.erase(u.begin() + 10000, u.end());
        v.clear();
    }

    if (w.empty()) {
        // found no successful path, so return an empty vector
        return npath;
    } else {
        BOOST_LOG_TRIVIAL(debug) << "have " << w.size() << " candidates";
        // find the most exact match, the one which covers all sequence with minimal
        // extra to end of graph, or longest
        auto longest_length = 0;
        std::vector<LocalNodePtr> longest_path;
        for (const auto& p : w) {
            candidate_string = "";
            for (const auto& s : p) {
                candidate_string += s->seq;
            }
            if (strcasecmp(query_string.c_str(), candidate_string.c_str()) == 0) {
                return p;
            } else if (candidate_string.size() > (uint)longest_length) {
                longest_path = p;
                longest_length = candidate_string.size();
            }
        }
        return longest_path;
    }
}

std::vector<LocalNodePtr> LocalGraph::top_path() const
{
    const bool graph_is_empty = nodes.empty();
    if (graph_is_empty) {
        fatal_error("Error getting top path in the graph: graph is empty");
    }

    std::vector<LocalNodePtr> npath;

    npath.push_back(nodes.at(0));
    while (not npath.back()->outNodes.empty()) {
        npath.push_back(npath.back()->outNodes[0]);
    }

    return npath;
}

std::vector<LocalNodePtr> LocalGraph::bottom_path() const
{
    const bool graph_is_empty = nodes.empty();
    if (graph_is_empty) {
        fatal_error("Error getting bottom path in the graph: graph is empty");
    }

    std::vector<LocalNodePtr> npath;

    npath.push_back(nodes.at(0));
    while (!npath.back()->outNodes.empty()) {
        npath.push_back(npath.back()->outNodes.back());
    }

    return npath;
}

bool LocalGraph::operator==(const LocalGraph& y) const
{
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {
        return false;
    }

    // false if have different nodes
    for (const auto& c : nodes) {
        // if node id doesn't exist
        auto it = y.nodes.find(c.first);
        if (it == y.nodes.end()) {
            return false;
        }
        // or node entries are different
        if (!(*c.second == *(it->second))) {
            return false;
        }
    }
    return true;
}

bool LocalGraph::operator!=(const LocalGraph& y) const { return !(*this == y); }

std::ostream& operator<<(std::ostream& out, LocalGraph const& data)
{
    for (const auto& c : data.nodes) {
        out << c.second->id << std::endl;
        for (uint32_t j = 0; j != c.second->outNodes.size(); ++j) {
            out << c.second->id << "->" << c.second->outNodes[j]->id << std::endl;
        }
    }
    return out;
}
