#include <stdint.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <cstring>
#include <cmath>
#include <utility>
#include "inthash.h"

using namespace std;

/* Taken from Heng Li minimap https://github.com/lh3/minimap/blob/master/sketch.c
 *
 * Licence for this repository:
 *
 * The MIT License
 *
 * Copyright (c) 2015 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * Note that this uses Thomas Wang's integer hash functions. 
 * See <https://gist.github.com/lh3/59882d6b96166dfc3d8d> for a snapshot.
 * Also for any 1<k<=64, let mask=(1<<k)-1. hash_64() is a bijection on [0,1<<k), which means
 *  hash_64(x, mask)==hash_64(y, mask) if and only if x==y.
 */

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void test_table()
{
    char c = 'A';
    int u = seq_nt4_table[(uint8_t)c];
    assert(u==0);
    c = 'a';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==0);

    c = 'C';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==1);
    c = 'c';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==1);

    c = 'G';
    u = seq_nt4_table[(uint8_t)c];
    assert(u=2);
    c = 'g';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==2);

    c = 'T';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==3);
    c = 't';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==3);

    c = 'N';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    c = 'R';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    c = 'Y';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    c = 'X';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    c = 'S';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    c = '-';
    u = seq_nt4_table[(uint8_t)c];
    assert(u==4);
    return;
}

uint64_t hash64(uint64_t& key, const uint64_t& mask)
{
	assert(key<=mask);
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

/* Now use these functions in my own code */

pair<uint64_t, uint64_t> kmerhash(const std::string& s, const uint32_t k)
{
    // this takes the hash of both forwards and reverse complement kmers and returns them as a pair 
    assert(s.size() == k);
    int c;
    uint64_t mask = pow(4,k) - 1, kh = 0, rckh = 0;
    char myArray[s.size()+1];//as 1 char space for null is also required
    strcpy(myArray, s.c_str());
    for (uint32_t i = 0; i < s.size(); ++i)
    {
        c = seq_nt4_table[(uint8_t)s[i]];
        //cout << s[i] << " -> " << c;
        if (c < 4) { // not an ambiguous base
	    kh += c*pow(4,i);
	    rckh += (3-c)*pow(4,k-i-1);
	}
        //cout << " so kh is now " << kh << endl;
    } 
    //cout << "s: " << s << " -> " << kh << endl;
    kh = hash64(kh, mask);
    rckh = hash64(rckh, mask);
    return make_pair(kh, rckh);
}


