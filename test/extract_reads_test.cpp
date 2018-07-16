#include <iostream>
#include <memory>
#include <set>
#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "extract_reads.h"
#include "interval.h"
#include "localPRG.h"
#include "minihit.h"
#include "minihits.h"
#include "prg/path.h"
#include "pangenome/pannode.h"
#include "pangenome/panread.h"

typedef std::shared_ptr<pangenome::Node> PanNodePtr;
typedef std::shared_ptr<pangenome::Read> PanReadPtr;

using namespace std;

TEST(ExtractReadsTest, identify_regions_null) {
    vector<uint32_t> covgs;
    auto low_covg_intervals = identify_regions(covgs, 0);
    vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2);
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_1) {
    vector<uint32_t> covgs({5, 6, 7, 5, 6, 6, 4, 4, 3, 5, 1, 1, 2, 3, 2, 4, 3});
    auto low_covg_intervals = identify_regions(covgs, 0, 0);
    vector<Interval> expected_intervals;
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 0);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 1);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 2);
    expected_intervals = {Interval(10, 12)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 1, 3);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_2) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 2, 0);
    vector<Interval> expected_intervals = {Interval(10,13),Interval(14,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 1);
    expected_intervals = {Interval(10,13),Interval(14,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 2);
    expected_intervals = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 3);
    expected_intervals = {Interval(10,13)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 2, 4);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_3) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 3, 0);
    vector<Interval> expected_intervals = {Interval(8,9),Interval(10,15), Interval(16,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 1);
    expected_intervals = {Interval(8,9),Interval(10,15), Interval(16,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 2);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 3);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 4);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 5);
    expected_intervals = {Interval(10,15)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 3, 6);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, identify_regions_4) {
    vector<uint32_t> covgs({5,6,7,5,6,6,4,4,3,5,1,1,2,3,2,4,3});
    auto low_covg_intervals = identify_regions(covgs, 4, 0);
    vector<Interval> expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 1);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 2);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 3);
    expected_intervals = {Interval(6,9),Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 4);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 5);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 6);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 7);
    expected_intervals = {Interval(10,17)};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);

    low_covg_intervals = identify_regions(covgs, 4, 8);
    expected_intervals = {};
    EXPECT_ITERABLE_EQ(vector<Interval>, expected_intervals, low_covg_intervals);
}

TEST(ExtractReadsTest, find_interval_in_localpath_minimal) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,3), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[2]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_short) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(1,4), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_uneven) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,4), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);

    found_path = find_interval_in_localpath(Interval(1,3), lmp);

    exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest, find_interval_in_localpath_multiple_sites) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {l3.prg.nodes[0], l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7], l3.prg.nodes[9]};
    // A G C T CGG  TAT
    auto found_path = find_interval_in_localpath(Interval(2,5), lmp);

    vector<LocalNodePtr> exp_path = {l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                     l3.prg.nodes[6], l3.prg.nodes[7]};
    EXPECT_ITERABLE_EQ(vector<LocalNodePtr>, exp_path, found_path);
}

TEST(ExtractReadsTest,hits_along_path) {
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {//l3.prg.nodes[0],
                                l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
                                };
    // A G C T CGG  TAT
    for (auto n : lmp){
        cout << n->pos << " ";
    }
    cout << endl;

    set<MinimizerHitPtr, pComp_path> hits, expected_subset;
    uint32_t read_id = 0, prg_id = 3, knode_id = 0;
    Interval read_interval(1,4);
    bool orientation(true);
    deque<Interval> d = {Interval(7,8), Interval(10, 12)};
    Path prg_path;
    prg_path.initialize(d);

    // hit not on path
    MinimizerHitPtr mh(make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation));
    hits.insert(mh);

    // hits branching from path
    d = {Interval(7,8), Interval(16, 17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31,33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28,30), Interval(33,33), Interval(40,41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4,5), Interval(8,9), Interval(16,17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(8,9), Interval(16,17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(16,17), Interval(27,29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);
    d = {Interval(27,30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, read_interval, prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    expected_subset.insert(mh);

    // hit against different prg
    // hit different orientation
    set<MinimizerHitPtr, pComp_path> found_subset = hits_along_path(hits,lmp);
    EXPECT_EQ(expected_subset.size(), found_subset.size());
    auto jt=found_subset.begin();
    for (auto it=expected_subset.begin(); it!=expected_subset.end() and jt!=found_subset.end();++it)
    {
        EXPECT_EQ(**jt, **it);
        ++jt;
    }

}

TEST(ExtractReadsTest, get_read_overlap_coordinates) {
    //
    //  Read 0 has prg 3 sequence in interval (2,12] only
    //  Read 1 has prg 3 sequence in interval (6,16] as well as noise
    //  Read 2 has prg 3 sequence in interval (4,20] stretched out
    //  Read 3 has prg 3 sequence in interval (4,14] but is missing bits
    //  Read 4 doesn't have prg 3 sequence - on all hits are noise
    //
    uint32_t pnode_id=3, prg_id=3, read_id = 0, knode_id = 0;
    string pnode_name = "three";
    bool orientation(true);
    deque<Interval> d;
    Path prg_path;
    MinimizerHitPtr mh;

    PanNodePtr pn= make_shared<pangenome::Node>(pnode_id, prg_id, pnode_name);
    PanReadPtr pr = make_shared<pangenome::Read>(read_id);

    set<MinimizerHitPtr, pComp> hits;


    // READ 0
    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(2,5), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28,30), Interval(33,33), Interval(40,41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7,10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4,5), Interval(8,9), Interval(16,17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(3,6), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8,9), Interval(16,17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4,7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16,17), Interval(27,29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5,8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27,30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6,9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id,hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 1
    read_id = 1;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6,9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(12,15), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28,30), Interval(33,33), Interval(40,41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(11,14), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4,5), Interval(8,9), Interval(16,17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7,10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8,9), Interval(16,17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16,17), Interval(27,29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27,30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10,13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7,8), Interval(16, 17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1,4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31,33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78,81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13,16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id,hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 2
    read_id = 2;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4,7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17,20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28,30), Interval(33,33), Interval(40,41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(15,18), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(4,5), Interval(8,9), Interval(16,17)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(5,8), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(8,9), Interval(16,17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16,17), Interval(27,29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(27,30)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10,13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7,8), Interval(16, 17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1,4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31,33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78,81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13,16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id,hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 3
    read_id = 3;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4,7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(10,13), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(28,30), Interval(33,33), Interval(40,41)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // hits on path
    d = {Interval(8,9), Interval(16,17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(6,9), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(16,17), Interval(27,29)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(7,10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);


    // noise
    d = {Interval(7,8), Interval(16, 17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1,4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(7,10), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id,hits);
    pn->reads.insert(pr);
    hits.clear();

    // READ 4
    read_id = 4;
    pr = make_shared<pangenome::Read>(read_id);

    // hits overlapping edges of path
    d = {Interval(0,1), Interval(4,5), Interval(8,9)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(4,7), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29,30), Interval(33,33), Interval(40,42)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(17,20), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    // noise
    d = {Interval(7,8), Interval(16, 17), Interval(27,28)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(1,4), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    mh = make_shared<MinimizerHit>(read_id, Interval(8,11), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(29, 30), Interval(31,33)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(9,12), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);
    d = {Interval(78,81)};
    prg_path.initialize(d);
    mh = make_shared<MinimizerHit>(read_id, Interval(13,16), prg_id, prg_path, knode_id, orientation);
    hits.insert(mh);

    pr->add_hits(prg_id,hits);
    pn->reads.insert(pr);
    hits.clear();

    // RUN GET_READ_OVERLAPS
    LocalPRG l3(3,"nested varsite", "A 5 G 7 C 8 T 7 T 9 CCG 10 CGG 9  6 G 5 TAT");
    vector<LocalNodePtr> lmp = {//l3.prg.nodes[0],
                                l3.prg.nodes[1], l3.prg.nodes[2], l3.prg.nodes[4],
                                l3.prg.nodes[6], l3.prg.nodes[7]//, l3.prg.nodes[9]
                                };
    // A G C T CGG  TAT

    vector<vector<uint32_t>> overlaps;
    vector<vector<uint32_t>> expected_overlaps = {{0, 3, 9, 1}, {1, 7, 13, 1}, {2, 5, 13, 1}, {3, 6, 10, 1}};

    get_read_overlap_coordinates(pn, overlaps, lmp);

    EXPECT_EQ(expected_overlaps.size(), overlaps.size());

    uint j = 0;
    for (auto coord : overlaps)
    {
        EXPECT_ITERABLE_EQ(vector<uint32_t>, expected_overlaps[j], coord);
        j++;
    }
}
