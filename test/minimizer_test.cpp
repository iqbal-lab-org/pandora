//#include <limits.h>
#include "gtest/gtest.h"
#include "minimizer.h"
#include "path.h"
#include "interval.h"
#include <set>
#include <vector>
#include <stdint.h> 

using std::set;
using namespace std;

class MinimizerTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MinimizerTest,comparisonCheck){
    list<interval> v1 = {interval(0,15)};
    list<interval> v2 = {interval(1,12), interval(14,18)};
    list<interval> v3 = {interval(1,16)};
    list<interval> v4 = {interval(0,10), interval(15,20)};

    Minimizer m1 = Minimizer("abcde", v1);
    Minimizer m2 = Minimizer("abcdg", v1);
    Minimizer m3 = Minimizer("abcde", v2);
    Minimizer m4 = Minimizer("abcde", v3);
    Minimizer m5 = Minimizer("abcde", v4);
    
    set<Minimizer> s;
    s.insert(m1);
    s.insert(m2);
    s.insert(m3);
    s.insert(m4);
    s.insert(m5);

    uint32_t j = 5;
    EXPECT_EQ(s.size(),j);

    vector<Minimizer> v = {m1, m5, m4, m3, m2};
    int i = 0;
    for (std::set<Minimizer>::iterator it=s.begin(); it!=s.end(); ++it)
    {
	it->print();
	v[i].print();
	EXPECT_EQ(it->miniWord, v[i].miniWord);
	//EXPECT_EQ(it->startPosOnString, v[i].startPosOnString);
	//EXPECT_EQ(it->endPosOnString, v[i].endPosOnString);
    	++i;
    }
}

int main(int ac, char* av[])
{
  testing::InitGoogleTest(&ac, av);
  return RUN_ALL_TESTS();
}
