#include <iostream>
#include "interval.h"
#include "gene_interval_info.h"


bool GeneIntervalInfo::operator<(const GeneIntervalInfo &y) const {
    cout << "check if pnodes defined" << endl;
    if (!pnode and !y.pnode) { return false; }

    cout << "check if one is defined" << endl;
    if (pnode and !y.pnode) { return false; }
    cout << "again" << endl;
    if (y.pnode and !pnode) {
        cout << "pnode not defined" << endl;
        return true;
    }

    cout << "continue checking" << endl;
    if (pnode->prg_id < y.pnode->prg_id) { return true; }
    if (y.pnode->prg_id < pnode->prg_id) { return false; }

    if (interval < y.interval) { return true; }
    if (y.interval < interval) { return false; }

    if (pnode->node_id < y.pnode->node_id) { return true; }
    if (y.pnode->node_id < pnode->node_id) { return false; }

    return false;
}

bool GeneIntervalInfo::operator==(const GeneIntervalInfo &y) const {
    if (!(interval == y.interval)) { return false; }
    if (!pnode and !y.pnode) { return true; }
    if (!pnode or !y.pnode) { return false; }
    if (pnode->prg_id != y.pnode->prg_id) { return false; }
    if (!(interval == y.interval)) { return false; }
    if (pnode->node_id != y.pnode->node_id) { return false; }

    return true;
}

bool GeneIntervalInfo::operator!=(const GeneIntervalInfo &y) const {
    return !(*this == y);
}

std::ostream &operator<<(std::ostream &out, GeneIntervalInfo const &y) {
    out << "{" << y.pnode->prg_id << ", " << y.interval << ", " << y.seq << "}";
    return out;
}
