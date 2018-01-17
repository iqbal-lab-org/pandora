#include "gtest/gtest.h"
#include "test_macro.cpp"
#include "noise_filtering.h"
#include "pangenome_graph_class.h"
#include "minihit.h"

using namespace std;

class NoiseFilteringTest : public ::testing::Test {
protected:
    virtual void SetUp() {
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
    }
};

TEST_F(NoiseFilteringTest, node_plus_orientation_to_num){
    EXPECT_EQ((uint16_t)0, node_plus_orientation_to_num(0,false));
    EXPECT_EQ((uint16_t)1, node_plus_orientation_to_num(0,true));
    EXPECT_EQ((uint16_t)2, node_plus_orientation_to_num(1,false));
    EXPECT_EQ((uint16_t)3, node_plus_orientation_to_num(1,true));
}

TEST_F(NoiseFilteringTest, num_to_node_plus_orientation){
    uint16_t node_id;
    bool node_orientation;

    num_to_node_plus_orientation(node_id, node_orientation, 0);
    EXPECT_EQ((uint16_t)0, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 1);
    EXPECT_EQ((uint16_t)0, node_id);
    EXPECT_EQ(node_orientation, true);

    num_to_node_plus_orientation(node_id, node_orientation, 2);
    EXPECT_EQ((uint16_t)1, node_id);
    EXPECT_EQ(node_orientation, false);

    num_to_node_plus_orientation(node_id, node_orientation, 3);
    EXPECT_EQ((uint16_t)1, node_id);
    EXPECT_EQ(node_orientation, true);
}

TEST_F(NoiseFilteringTest, rc_num)
{
    EXPECT_EQ(node_plus_orientation_to_num(0,false), rc_num(node_plus_orientation_to_num(0,true)));
    EXPECT_EQ(node_plus_orientation_to_num(0,true), rc_num(node_plus_orientation_to_num(0,false)));
    EXPECT_EQ(node_plus_orientation_to_num(1,false), rc_num(node_plus_orientation_to_num(1,true)));
    EXPECT_EQ(node_plus_orientation_to_num(1,true), rc_num(node_plus_orientation_to_num(1,false)));
    EXPECT_EQ(node_plus_orientation_to_num(2,false), rc_num(node_plus_orientation_to_num(2,true)));
    EXPECT_EQ(node_plus_orientation_to_num(2,true), rc_num(node_plus_orientation_to_num(2,false)));
}

TEST_F(NoiseFilteringTest,hashed_node_ids_to_ids_and_orientations)
{
    deque<uint16_t> d = {3,1,2,0};
    vector<uint16_t> v;
    vector<uint16_t> v_exp = {1,0,1,0};
    vector<bool> b;
    vector<bool> b_exp = {true, true, false, false};

    hashed_node_ids_to_ids_and_orientations(d,v,b);
    EXPECT_ITERABLE_EQ(vector<uint16_t>, v_exp, v);
    EXPECT_ITERABLE_EQ(vector<bool>, b_exp, b);
}

TEST_F(NoiseFilteringTest,overlap_forwards)
{
    deque<uint16_t> d1 = {0,1,2};
    deque<uint16_t> d2 = {1,2,3};
    EXPECT_EQ(overlap_forwards(d1,d2), true);
    EXPECT_EQ(overlap_forwards(d2,d1), false);
    EXPECT_EQ(overlap_forwards(d1,d1), false);
    EXPECT_EQ(overlap_forwards(d2,d2), false);

    // works when d1 longer than d2
    d1 = {0,4,6,2,5,4,0,1,2};
    EXPECT_EQ(overlap_forwards(d1,d2), true);
    EXPECT_EQ(overlap_forwards(d1,d1), false);

    // if overlap > 1 then is false
    d2 = {1,2,3,4};
    EXPECT_EQ(overlap_forwards(d1,d2), false);
    EXPECT_EQ(overlap_forwards(d2,d2), false);

    // if d2 longer than d1, should be false
    d2 = {0,4,6,2,5,4,0,1,2,3};
    EXPECT_DEATH(overlap_forwards(d1,d2), "");
    EXPECT_EQ(overlap_forwards(d2,d1), false);
}

TEST_F(NoiseFilteringTest,overlap_backwards)
{
    deque<uint16_t> d1 = {0,1,2};
    deque<uint16_t> d2 = {1,2,3};
    EXPECT_EQ(overlap_backwards(d2,d1), true);
    EXPECT_EQ(overlap_backwards(d1,d2), false);
    EXPECT_EQ(overlap_backwards(d1,d1), false);
    EXPECT_EQ(overlap_backwards(d2,d2), false);

    // works when d1 longer than d2
    d1 = {0,4,6,2,5,4,0,1,2};
    d2 = {1,0,4};
    EXPECT_EQ(overlap_backwards(d1,d2), true);
    EXPECT_EQ(overlap_backwards(d1,d1), false);
    EXPECT_EQ(overlap_backwards(d2,d2), false);


    // if overlap > 1 then is false
    d2 = {1,2,0,4};
    EXPECT_EQ(overlap_backwards(d1,d2), false);
    EXPECT_EQ(overlap_backwards(d2,d2), false);

    // if d2 longer than d1, should still work?
    d2 = {3,0,4,6,2,5,4,0,1,2,4,6};
    EXPECT_EQ(overlap_backwards(d1,d2), true);
}

TEST_F(NoiseFilteringTest,reverse_hashed_node)
{
    deque<uint16_t> d1 = {0,1,2};
    deque<uint16_t> d2 = {3,0,1};
    EXPECT_ITERABLE_EQ(deque<uint16_t>, d1, reverse_hashed_node(d2));
    EXPECT_ITERABLE_EQ(deque<uint16_t>, d2, reverse_hashed_node(d1));
}

TEST_F(NoiseFilteringTest,dbg_node_ids_to_ids_and_orientations)
{
    deque<uint16_t> read1 = {0,1,3,5,6};
    deque<bool> read1_o = {0,0,0,1,0};
    deque<uint16_t> read2 = {3,7,6,5,3,1};
    deque<bool> read2_o = {0,0,1,0,1,1};

    vector<uint16_t> node_ids;
    vector<bool> node_orients;
    vector<uint16_t> exp = {0,1,3,5,6,7,3};
    vector<bool> exp_o = {0,0,0,1,0,1,1};

    debruijn::Graph dbg(3);
    deque<uint16_t> d = {0,2,6};
    dbg.add_node(d, 0);
    d = {2,6,11};
    dbg.add_node(d,0);
    d = {6,11,12};
    dbg.add_node(d,0);

    d = {6,14,13};
    dbg.add_node(d,1);
    d = {14,13,10};
    dbg.add_node(d,1);
    d = {13,10,7};
    dbg.add_node(d,1);
    d = {10,7,3};
    dbg.add_node(d,1);

    EXPECT_EQ(dbg.nodes.size(), (uint)5);
    deque<uint16_t> tig = {0,1,2,4,3};
    dbg_node_ids_to_ids_and_orientations(dbg, tig, node_ids, node_orients);
    EXPECT_ITERABLE_EQ(vector<uint16_t>, exp, node_ids);
    EXPECT_ITERABLE_EQ(vector<bool>, exp_o, node_orients);

    // also want to find if tig is listed in reverse

}

TEST_F(NoiseFilteringTest,construct_debruijn_graph_from_pangraph)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(4,"4",0, mhs);
    pg.add_node(5,"5",0, mhs);

    // overlaps to create loop
    pg.add_node(3,"3",1, mhs);
    pg.add_node(4,"4",1, mhs);
    pg.add_node(5,"5",1, mhs);
    pg.add_node(0,"0",1, mhs);
    pg.add_node(1,"1",1, mhs);
    pg.add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg.add_node(1,"1",2, mhs);
    pg.add_node(2,"2",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(7,"7",2, mhs);

    // all disjoint, short
    pg.add_node(0,"0",3, mhs);
    pg.add_node(6,"6",3, mhs);
    pg.add_node(3,"3",3, mhs);
    pg.add_node(4,"4",3, mhs);

    // deviates in middle
    pg.add_node(0,"0",4, mhs);
    pg.add_node(1,"1",4, mhs);
    pg.add_node(2,"2",4, mhs);
    pg.add_node(6,"6",4, mhs);
    pg.add_node(3,"3",4, mhs);
    pg.add_node(4,"4",4, mhs);
    pg.add_node(5,"5",4, mhs);

    // all disjoint, long
    pg.add_node(6,"6",5, mhs);
    pg.add_node(1,"1",5, mhs);
    pg.add_node(2,"2",5, mhs);
    pg.add_node(6,"6",5, mhs);
    pg.add_node(3,"3",5, mhs);

    debruijn::Graph dbg = construct_debruijn_graph_from_pangraph(3,pg);

    debruijn::Graph dbg_exp(3);
    deque<uint16_t> d = {0,2,4};
    debruijn::NodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2,4,6};
    debruijn::NodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {4,6,8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);

    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);
    d = {10,0,2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1,n2);
    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,6};
    n1 = dbg_exp.add_node(d, 2);
    d = {4,6,14};
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1,n2);

    d = {0,12,6};
    n1 = dbg_exp.add_node(d, 3);
    d = {12,6,8};
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1,n2);

    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 4);
    d = {2,4,12};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {4,12,6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);
    d = {12,6,8};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {6,8,10};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);

    d = {12,2,4};
    n1 = dbg_exp.add_node(d, 5);
    d = {2,4,12};
    n2 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n1,n2);
    d = {4,12,6};
    n1 = dbg_exp.add_node(d, 5);
    dbg_exp.add_edge(n2,n1);

    EXPECT_EQ(dbg_exp, dbg);
}

TEST_F(NoiseFilteringTest,remove_leaves)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(4,"4",0, mhs);
    pg.add_node(5,"5",0, mhs);

    // overlapping in loop
    pg.add_node(3,"3",1, mhs);
    pg.add_node(4,"4",1, mhs);
    pg.add_node(5,"5",1, mhs);
    pg.add_node(0,"0",1, mhs);
    pg.add_node(1,"1",1, mhs);
    pg.add_node(2,"2",1, mhs);


    // starts correct and deviates
    pg.add_node(1,"1",2, mhs);
    pg.add_node(2,"2",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(7,"7",2, mhs);

    // incorrect short
    pg.add_node(0,"0",3, mhs);
    pg.add_node(5,"5",3, mhs);//6
    pg.add_node(3,"3",3, mhs);
    pg.add_node(4,"4",3, mhs);

    // deviates in middle
    pg.add_node(0,"0",4, mhs);
    pg.add_node(1,"1",4, mhs);
    pg.add_node(2,"2",4, mhs);
    pg.add_node(6,"6",4, mhs);
    pg.add_node(3,"3",4, mhs);
    pg.add_node(4,"4",4, mhs);
    pg.add_node(5,"5",4, mhs);

    // incorrect longer
    pg.add_node(6,"6",5, mhs);
    pg.add_node(1,"1",5, mhs);
    pg.add_node(1,"1",5, mhs);//2
    pg.add_node(6,"6",5, mhs);
    pg.add_node(3,"3",5, mhs);

    cout << "pg is now: " << endl << pg << endl;

    debruijn::Graph dbg = construct_debruijn_graph_from_pangraph(3,pg);
    remove_leaves(pg, dbg);

    debruijn::Graph dbg_exp(3);
    deque<uint16_t> d = {0,2,4};
    debruijn::NodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2,4,6};
    debruijn::NodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {4,6,8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);

    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);
    d = {10,0,2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1,n2);
    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,6};
    n1 = dbg_exp.add_node(d, 2);

    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 4);
    d = {2,4,12};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {4,12,6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);
    d = {12,6,8};
    n2 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n1,n2);
    d = {6,8,10};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,12};
    n2 = dbg_exp.add_node(d, 4);
    d = {4,12,6};
    n1 = dbg_exp.add_node(d, 4);
    dbg_exp.add_edge(n2,n1);

    EXPECT_EQ(dbg_exp, dbg);

    pangenome::Graph pg_exp;
    pg_exp.add_node(0,"0",0, mhs);
    pg_exp.add_node(1,"1",0, mhs);
    pg_exp.add_node(2,"2",0, mhs);
    pg_exp.add_node(3,"3",0, mhs);
    pg_exp.add_node(4,"4",0, mhs);
    pg_exp.add_node(5,"5",0, mhs);

    pg_exp.add_node(3,"3",1,mhs);
    pg_exp.add_node(4,"4",1, mhs);
    pg_exp.add_node(5,"5",1, mhs);
    pg_exp.add_node(0,"0",1, mhs);
    pg_exp.add_node(1,"1",1, mhs);
    pg_exp.add_node(2,"2",1, mhs);

    pg_exp.add_node(1,"1",2, mhs);
    pg_exp.add_node(2,"2",2, mhs);
    pg_exp.add_node(3,"3",2, mhs);

    pg_exp.add_node(0,"0",4, mhs);
    pg_exp.add_node(1,"1",4, mhs);
    pg_exp.add_node(2,"2",4, mhs);
    pg_exp.add_node(6,"6",4, mhs);
    pg_exp.add_node(3,"3",4, mhs);
    pg_exp.add_node(4,"4",4, mhs);
    pg_exp.add_node(5,"5",4, mhs);

    EXPECT_EQ(pg_exp, pg);
}

TEST_F(NoiseFilteringTest,filter_unitigs)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(4,"4",0, mhs);
    pg.add_node(5,"5",0, mhs);
    pg.add_node(0,"0",0, mhs);

    // overlapping in loop
    pg.add_node(3,"3",1, mhs);
    pg.add_node(4,"4",1, mhs);
    pg.add_node(5,"5",1, mhs);
    pg.add_node(0,"0",1, mhs);
    pg.add_node(1,"1",1, mhs);
    pg.add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg.add_node(1,"1",2, mhs);
    pg.add_node(2,"2",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(7,"7",2, mhs);

    // incorrect short
    pg.add_node(0,"0",3, mhs);
    pg.add_node(5,"5",3, mhs);//6
    pg.add_node(3,"3",3, mhs);
    pg.add_node(4,"4",3, mhs);

    // deviates in middle
    pg.add_node(0,"0",4, mhs);
    pg.add_node(1,"1",4, mhs);
    pg.add_node(2,"2",4, mhs);
    pg.add_node(6,"6",4, mhs);
    pg.add_node(3,"3",4, mhs);
    pg.add_node(4,"4",4, mhs);
    pg.add_node(5,"5",4, mhs);

    cout << "original pg is: " << endl << pg << endl;

    debruijn::Graph dbg = construct_debruijn_graph_from_pangraph(3,pg);
    filter_unitigs(pg, dbg, 1);

    cout << "dbg is now: " << endl << dbg << endl;

    debruijn::Graph dbg_exp(3);
    deque<uint16_t> d = {0,2,4};
    debruijn::NodePtr n1 = dbg_exp.add_node(d, 0);
    d = {2,4,6};
    debruijn::NodePtr n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {4,6,8};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n1,n2);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 0);
    dbg_exp.add_edge(n2,n1);

    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 1);
    d = {8,10,0};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);
    d = {10,0,2};
    n2 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n1,n2);
    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 1);
    dbg_exp.add_edge(n2,n1);

    d = {2,4,6};
    n1 = dbg_exp.add_node(d, 2);
    d = {4,6,14};
    n2 = dbg_exp.add_node(d, 2);
    dbg_exp.add_edge(n1,n2);

    d = {0,10,6};
    n1 = dbg_exp.add_node(d, 3);
    d = {10,6,8};
    n2 = dbg_exp.add_node(d, 3);
    dbg_exp.add_edge(n1,n2);

    d = {0,2,4};
    n1 = dbg_exp.add_node(d, 4);
    d = {6,8,10};
    n2 = dbg_exp.add_node(d, 4);

    EXPECT_EQ(dbg_exp, dbg);

    cout << "pg is now: " << endl << pg << endl;

    pangenome::Graph pg_exp;
    pg_exp.add_node(0,"0",0, mhs);
    pg_exp.add_node(1,"1",0, mhs);
    pg_exp.add_node(2,"2",0, mhs);
    pg_exp.add_node(3,"3",0, mhs);
    pg_exp.add_node(4,"4",0, mhs);
    pg_exp.add_node(5,"5",0, mhs);
    pg_exp.add_node(0,"0",0, mhs);

    pg_exp.add_node(3,"3",1,mhs);
    pg_exp.add_node(4,"4",1, mhs);
    pg_exp.add_node(5,"5",1, mhs);
    pg_exp.add_node(0,"0",1, mhs);
    pg_exp.add_node(1,"1",1, mhs);
    pg_exp.add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg_exp.add_node(1,"1",2, mhs);
    pg_exp.add_node(2,"2",2, mhs);
    pg_exp.add_node(3,"3",2, mhs);
    pg_exp.add_node(7,"7",2, mhs);

    // incorrect short
    pg_exp.add_node(0,"0",3, mhs);
    pg_exp.add_node(5,"5",3, mhs);//6
    pg_exp.add_node(3,"3",3, mhs);
    pg_exp.add_node(4,"4",3, mhs);

    // deviates in middle
    pg_exp.add_node(0,"0",4, mhs);
    pg_exp.add_node(1,"1",4, mhs);
    pg_exp.add_node(2,"2",4, mhs);
    pg_exp.add_node(3,"3",4, mhs);
    pg_exp.add_node(4,"4",4, mhs);
    pg_exp.add_node(5,"5",4, mhs);

    cout << "so google test" << endl;
    EXPECT_EQ(pg_exp, pg);
}


TEST_F(NoiseFilteringTest,detangle_pangraph_with_debruijn_graph)
{
    set<MinimizerHitPtr, pComp> mhs;
    pangenome::Graph pg;
    pg.add_node(0,"0",0, mhs);
    pg.add_node(1,"1",0, mhs);
    pg.add_node(2,"2",0, mhs);
    pg.add_node(3,"3",0, mhs);
    pg.add_node(4,"4",0, mhs);
    pg.add_node(5,"5",0, mhs);
    pg.add_node(0,"0",0, mhs);

    // overlapping in loop
    pg.add_node(3,"3",1, mhs);
    pg.add_node(4,"4",1, mhs);
    pg.add_node(5,"5",1, mhs);
    pg.add_node(0,"0",1, mhs);
    pg.add_node(1,"1",1, mhs);
    pg.add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg.add_node(1,"1",2, mhs);
    pg.add_node(2,"2",2, mhs);
    pg.add_node(3,"3",2, mhs);
    pg.add_node(7,"7",2, mhs);

    // incorrect short
    pg.add_node(0,"0",3, mhs);
    pg.add_node(5,"5",3, mhs);//6
    pg.add_node(3,"3",3, mhs);
    pg.add_node(4,"4",3, mhs);

    // deviates in middle
    pg.add_node(0,"0",4, mhs);
    pg.add_node(1,"1",4, mhs);
    pg.add_node(2,"2",4, mhs);
    pg.add_node(6,"6",4, mhs);
    pg.add_node(3,"3",4, mhs);
    pg.add_node(4,"4",4, mhs);
    pg.add_node(5,"5",4, mhs);

    cout << "original pg is: " << endl << pg << endl;

    debruijn::Graph dbg = construct_debruijn_graph_from_pangraph(3,pg);
    detangle_pangraph_with_debruijn_graph(pg, dbg);

    /*pangenome::Graph pg_exp;
    pg_exp.add_node(8,"0",0, mhs);
    pg_exp.add_node(9,"1",0, mhs);
    pg_exp.add_node(10,"2",0, mhs);
    pg_exp.add_node(11,"3",0, mhs);
    pg_exp.add_node(12,"4",0, mhs);
    pg_exp.add_node(13,"5",0, mhs);
    pg_exp.add_node(14,"0",0, mhs);

    pg_exp.add_node(3,"3",1,mhs);
    pg_exp.add_node(4,"4",1, mhs);
    pg_exp.add_node(5,"5",1, mhs);
    pg_exp.add_node(0,"0",1, mhs);
    pg_exp.add_node(1,"1",1, mhs);
    pg_exp.add_node(2,"2",1, mhs);

    // starts correct and deviates
    pg_exp.add_node(1,"1",2, mhs);
    pg_exp.add_node(2,"2",2, mhs);
    pg_exp.add_node(3,"3",2, mhs);
    pg_exp.add_node(7,"7",2, mhs);

    // incorrect short
    pg_exp.add_node(0,"0",3, mhs);
    pg_exp.add_node(5,"5",3, mhs);//6
    pg_exp.add_node(3,"3",3, mhs);
    pg_exp.add_node(4,"4",3, mhs);

    // deviates in middle
    pg_exp.add_node(0,"0",4, mhs);
    pg_exp.add_node(1,"1",4, mhs);
    pg_exp.add_node(2,"2",4, mhs);
    pg_exp.add_node(3,"3",4, mhs);
    pg_exp.add_node(4,"4",4, mhs);
    pg_exp.add_node(5,"5",4, mhs);*/
}