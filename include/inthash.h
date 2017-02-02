#ifndef __INTHASH_H_INCLUDED__   // if inthash.h hasn't been included yet...
#define __INTHASH_H_INCLUDED__

#include <cstdint>

uint64_t hash64(uint64_t& key, const uint64_t& mask);
std::pair<uint64_t,uint64_t> kmerhash(const std::string& s, const uint32_t k);
void test_table();

#endif
