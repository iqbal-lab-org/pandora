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