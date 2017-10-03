#include <iostream>
#include <cstring>
#include <string.h>
#include <map>
#include <fstream>
#include <cassert>
#include <vector>
#include "utils.h"
#include "localnode.h"
#include "localgraph.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

using namespace std;

LocalGraph::LocalGraph() {
    // reserve space in index
    //index.reserve(10);
}

LocalGraph::~LocalGraph() {
    for (auto c: nodes) {
        delete c.second;
    }
}

void LocalGraph::add_node(const uint32_t &id, const string &seq, const Interval &pos) {
    assert(seq.length() == pos.length);
    map<uint32_t, LocalNode *>::iterator it = nodes.find(id);
    if (it == nodes.end()) {
        LocalNode *n;
        n = new LocalNode(seq, pos, id);
        nodes[id] = n;
        //cout << "Added node " << id << endl;
    } else {
        assert((it->second->seq == seq) && (it->second->pos == pos));
    }
    return;
}

void LocalGraph::add_edge(const uint32_t &from, const uint32_t &to) {
    map<uint32_t, LocalNode *>::iterator from_it = nodes.find(from);
    map<uint32_t, LocalNode *>::iterator to_it = nodes.find(to);
    assert((from_it != nodes.end()) && (to_it != nodes.end()));
    if ((from_it != nodes.end()) && (to_it != nodes.end())) {
        LocalNode *f = (nodes.find(from)->second);
        LocalNode *t = (nodes.find(to)->second);
        assert(f->pos.end <= t->pos.start || assert_msg(
                f->pos.end << ">" << t->pos.start << " so cannot add edge from node " << *f << " to node " << *t));
        f->outNodes.push_back(t);
        //cout << "Added edge (" << f->id << ", " << t->id << ")" << endl;
    }
    return;
}

/*void LocalGraph::add_varsite (const uint8_t level, const uint32_t pre_site_id, const uint32_t post_site_id)
{
    assert(pre_site_id <= post_site_id);
    while (level >= index.size())
    {
        vector<pair<uint32_t, uint32_t>> levelv;
	levelv.reserve(400);
	//levelv = {};
	index.insert(index.end(), 1, levelv);
    }
    index[level].push_back(make_pair(pre_site_id, post_site_id));
    return;
}*/

void LocalGraph::write_gfa(const string &filepath) const {
    ofstream handle;
    handle.open(filepath);
    handle << "H\tVN:Z:1.0\tbn:Z:--linear --singlearr" << endl;
    for (map<uint32_t, LocalNode *>::const_iterator it = nodes.begin(); it != nodes.end(); ++it) {
        handle << "S\t" << it->second->id << "\t";
        if (it->second->seq == "") {
            handle << "*";
        } else {
            handle << it->second->seq;
        }
        handle << "\tRC:i:" << it->second->covg << endl;
        for (uint32_t j = 0; j < it->second->outNodes.size(); ++j) {
            handle << "L\t" << it->second->id << "\t+\t" << it->second->outNodes[j]->id << "\t+\t0M" << endl;
        }
    }
    handle.close();
}

void LocalGraph::read_gfa(const string &filepath) {
    uint32_t id, from, to;
    string line;
    vector<string> split_line;
    uint32_t i = 0;

    ifstream myfile(filepath);
    if (myfile.is_open()) {
        while (getline(myfile, line).good()) {
            if (line[0] == 'S') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 3);
                if (split_line[2] == "*") {
                    split_line[2] = "";
                }
                id = stoi(split_line[1]);
                add_node(id, (string) split_line[2], Interval(i, i + split_line[2].size()));
                i += split_line[2].size();
            }
        }

        myfile.clear();
        myfile.seekg(0, myfile.beg);
        while (getline(myfile, line).good()) {
            if (line[0] == 'L') {
                split_line = split(line, "\t");
                assert(split_line.size() >= 5);
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
        cerr << "Unable to open GFA file " << filepath << endl;
        exit(1);
    }
    return;
}

vector<Path> LocalGraph::walk(const uint32_t &node_id, const uint32_t &pos, const uint32_t &len) const {
    //cout << "walking graph from node " << node_id << " pos " << pos << " for length " << len << endl;
    // walks from position pos in node node for length len bases
    assert((nodes.at(node_id)->pos.start <= pos && nodes.at(node_id)->pos.end >= pos) || assert_msg(
            nodes.at(node_id)->pos.start << "<=" << pos << " and " << nodes.at(node_id)->pos.end << ">="
                                         << pos)); // if this fails, pos given lies on a different node
    vector<Path> return_paths, walk_paths;
    return_paths.reserve(20);
    walk_paths.reserve(20);
    Path p, p2;
    deque<Interval> d;

    //cout << "pos+len: " << pos+len << " nodes.at(node_id)->pos.end: " << nodes.at(node_id)->pos.end << endl;
    if (pos + len <= nodes.at(node_id)->pos.end) {
        d = {Interval(pos, pos + len)};
        p.initialize(d);
        //cout << "return path: " << p << endl;
        return_paths.push_back(p);
        //cout << "return_paths size: " << return_paths.size() << endl; 
        return return_paths;
    }
    uint32_t len_added = min(nodes.at(node_id)->pos.end - pos, len);

    //cout << "len: " << len << " len_added: " << len_added << endl;
    if (len_added < len) {
        for (vector<LocalNode *>::iterator it = nodes.at(node_id)->outNodes.begin();
             it != nodes.at(node_id)->outNodes.end(); ++it) {
            //cout << "Following node: " << (*it)->id << " to add " << len-len_added << " more bases" << endl;
            walk_paths = walk((*it)->id, (*it)->pos.start, len - len_added);
            //cout << "walk paths size: " << walk_paths.size() << endl;
            for (vector<Path>::iterator it2 = walk_paths.begin(); it2 != walk_paths.end(); ++it2) {
                // Note, would have just added start interval to each item in walk_paths, but can't seem to force result of it2 to be non-const
                //cout << (*it2) << endl;
                p2.initialize((*it2).path);
                p2.add_start_interval(Interval(pos, nodes.at(node_id)->pos.end));
                //cout << "path: " << p2 << " p2.length: " << p2.length << endl;
                if (p2.length() == len) {
                    return_paths.push_back(p2);
                }
            }
        }
    }
    return return_paths;
}

vector<Path> LocalGraph::walk_back(const uint32_t &node_id, const uint32_t &pos, const uint32_t &len) const {
    //cout << "start walking back from " << pos << " in node " << node_id << " for length " << len << endl;
    // walks from position pos in node back through prg for length len bases
    assert((nodes.at(node_id)->pos.start <= pos && nodes.at(node_id)->pos.end >= pos) || assert_msg(
            nodes.at(node_id)->pos.start << "<=" << pos << " and " << nodes.at(node_id)->pos.end << ">="
                                         << pos)); // if this fails, pos given lies on a different node
    vector<Path> return_paths, walk_paths;
    return_paths.reserve(20);
    walk_paths.reserve(20);
    Path p, p2;
    deque<Interval> d;

    if (nodes.at(node_id)->pos.start + len <= pos) {
        d = {Interval(pos - len, pos)};
        p.initialize(d);
        //cout << "return path: " << p << endl;
        return_paths.push_back(p);
        return return_paths;
    }

    uint32_t len_added = min(pos - nodes.at(node_id)->pos.start, len);
    //cout << "len: " << len << " len_added: " << len_added << endl;

    vector<LocalNode *>::iterator innode;
    if (len_added < len) {
        for (map<uint32_t, LocalNode *>::const_iterator it = nodes.begin(); it != nodes.find(node_id); ++it) {
            innode = find(it->second->outNodes.begin(), it->second->outNodes.end(), nodes.at(node_id));
            if (innode != it->second->outNodes.end()) {
                walk_paths = walk_back(it->second->id, it->second->pos.end, len - len_added);
                for (uint i = 0; i != walk_paths.size(); ++i) {
                    p2.initialize(walk_paths[i].path);
                    p2.add_end_interval(Interval(nodes.at(node_id)->pos.start, pos));
                    //cout << p2 << endl;
                    if (p2.length() == len) {
                        //cout << "output path: " << p2 << endl;
                        return_paths.push_back(p2);
                    }
                }
            }
        }
    }
    return return_paths;
}

LocalNode *LocalGraph::get_previous_node(const LocalNode *n) const {
    if (n->id == 0) {
        return nullptr;
    } else {
        for (auto c : nodes) {
            if (find(c.second->outNodes.begin(), c.second->outNodes.end(), n) != c.second->outNodes.end()) {
                return c.second;
            } else if (c.first > n->id) {
                break;
            }
        }
        // if we get here, there was no previous node to be found.
        return nullptr;
    }
}

vector<LocalNode *> LocalGraph::nodes_along_string(const string &query_string) const {
    // Note expects the query string to start at the start of the PRG - can change this later
    vector<vector<LocalNode *>> u, v;   // u <=> v
    // ie reject paths in u, or extend and add to v
    // then set u=v and continue
    u.reserve(100);
    v.reserve(100);
    vector<LocalNode *> npath;
    string candidate_string = "";
    bool extended = true;

    assert(nodes.size() > 0); //otherwise empty nodes -> segfault

    // if there is only one node in PRG, simple case, do simple string compare
    if (nodes.size() == 1 and strcasecmp(query_string.c_str(), nodes.at(0)->seq.c_str()) == 0) {
        return {nodes.at(0)};
    }

    u = {{nodes.at(0)}};

    while (u.size() > 0) {
        for (uint32_t i = 0; i != u.size(); ++i) {
            for (uint32_t j = 0; j != u[i].size(); ++j) {
                candidate_string += u[i][j]->seq;
            }

            for (uint32_t j = 0; j != u[i].back()->outNodes.size(); ++j) {
                // if the start of query_string matches extended candidate_string, want to query candidate path extensions
                //if ( query_string.substr(0,candidate_string.size()+u[i].back()->outNodes[j]->seq.size()) == candidate_string+u[i].back()->outNodes[j]->seq)
                if (strcasecmp(
                        query_string.substr(0, candidate_string.size() + u[i].back()->outNodes[j]->seq.size()).c_str(),
                        (candidate_string + u[i].back()->outNodes[j]->seq).c_str()) == 0) {
                    if (candidate_string.size() + u[i].back()->outNodes[j]->seq.size() >= query_string.size()) {
                        // we have now found the whole of the query_string
                        u[i].push_back(u[i].back()->outNodes[j]);
                        while (u[i].back()->outNodes.size() > 0 and extended == true) {
                            extended = false;
                            for (uint n = 0; n != u[i].back()->outNodes.size(); ++n) {
                                if (u[i].back()->outNodes[n]->pos.length == 0) {
                                    u[i].push_back(u[i].back()->outNodes[n]);
                                    extended = true;
                                    break;
                                }
                            }
                        }
                        return u[i];
                    } else {
                        v.push_back(u[i]);
                        v.back().push_back(u[i].back()->outNodes[j]);
                    }
                }
            }
            candidate_string = "";
        }
        u = v;
        v.clear();
    }
    // found no successful path, so return an empty vector
    return npath;
}

vector<LocalNode *> LocalGraph::top_path() const {
    vector<LocalNode *> npath;

    assert(nodes.size() > 0); //otherwise empty nodes -> segfault

    npath.push_back(nodes.at(0));
    while (npath.back()->outNodes.size() > 0) {
        npath.push_back(npath.back()->outNodes[0]);
    }

    return npath;
}

vector<LocalNode *> LocalGraph::bottom_path() const {
    vector<LocalNode *> npath;

    assert(nodes.size() > 0); //otherwise empty nodes -> segfault

    npath.push_back(nodes.at(0));
    while (npath.back()->outNodes.size() > 0) {
        npath.push_back(npath.back()->outNodes.back());
    }

    return npath;
}

bool LocalGraph::operator==(const LocalGraph &y) const {
    // false if have different numbers of nodes
    if (y.nodes.size() != nodes.size()) {//cout << "different numbers of nodes" << endl; 
        return false;
    }

    // false if have different nodes
    for (const auto c: nodes) {
        // if node id doesn't exist 
        map<uint32_t, LocalNode *>::const_iterator it = y.nodes.find(c.first);
        if (it == y.nodes.end()) {//cout << "node id doesn't exist" << endl;
            return false;
        }
        // or node entries are different
        if (!(*c.second ==
              *(it->second))) {//cout << "node id " << c.first << " exists but has different values" << endl;
            return false;
        }
    }
    // otherwise is true
    return true;
}

bool LocalGraph::operator!=(const LocalGraph &y) const {
    return !(*this == y);
}

std::ostream &operator<<(std::ostream &out, LocalGraph const &data) {
    for (const auto c: data.nodes) {
        out << c.second->id << endl;
        for (uint j = 0; j != c.second->outNodes.size(); ++j) {
            out << c.second->id << "->" << c.second->outNodes[j]->id << endl;
        }
    }
    return out;
}
