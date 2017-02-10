#include "gtest/gtest.h"
#include "inthash.h"
#include <iostream>
#include <cstring>
#include <string>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdint.h>
#include <iostream>

using namespace std;

class InthashTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(InthashTest,checkCharToInt){
    test_table();
}

set<string> generate_kmers(vector<string> v, uint32_t k)
{
    set<string> current_strings, new_strings;
    if (k >= 1)
    {
        for (uint32_t j=0;j<4;j++)
        {
            current_strings.insert(v[j]);
	}
    } else {
        return current_strings;
    }
    
    uint32_t n = 1;
    while ( n < k )
    {
	for (set<string>::iterator it = current_strings.begin(); it != current_strings.end(); ++it)
        {
   	    for (uint32_t j=0;j<4;j++)
            {
                new_strings.insert(*it + v[j]);
            }
        }
        current_strings = new_strings;
        new_strings.clear();
	n++;
    }
    return current_strings;
}

TEST_F(InthashTest,check1to1){
    vector<uint32_t> ks = {3,5};
    KmerHash hash;
    for (vector<uint32_t>::iterator jt = ks.begin(); jt != ks.end(); ++jt)
    {
        uint32_t k = *jt;
        // generate kmers for a given k
        vector<string> dna = {"A", "G", "T", "C"};
        set<string> kmers = generate_kmers(dna, k);
        // generate inthash of each kmer and assert it is different to the previous ones
        vector<uint64_t> khs, khs1, khs2;
        pair<uint64_t,uint64_t> kh;
        for (set<string>::iterator it = kmers.begin(); it != kmers.end(); ++it)
        {
            kh = hash.kmerhash(*it, k);
            EXPECT_EQ((kh.first < pow(4,k)), true);
            if (find(khs.begin(), khs.end(), kh.first) != khs.end())
            {
		cout << *it << ": " << kh.first << " == " << *find(khs.begin(), khs.end(), kh.first) << endl;
	    }
            EXPECT_EQ((find(khs.begin(), khs.end(), kh.first) == khs.end()), true);
            khs.push_back(kh.first); // so each kmerhash value is first in one

            // don't expect any kmerhash value to be in more than 2 pairs
	    EXPECT_EQ((find(khs2.begin(), khs2.end(), min(kh.first, kh.second)) == khs2.end()), true);
            if (find(khs1.begin(), khs1.end(), min(kh.first, kh.second)) != khs1.end())
	    {
		khs2.push_back(min(kh.first,kh.second)); // if we found it as min once before, add to things seen twice
	    } else {
		khs1.push_back(min(kh.first,kh.second)); // otherwise, add to list seen once
	    }
        }
    }
}

