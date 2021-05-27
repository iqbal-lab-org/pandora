#ifndef __INTHASH_H_INCLUDED__ // if inthash.h hasn't been included yet...
#define __INTHASH_H_INCLUDED__


#include <cstdint>
#include <string> //cstring doesn't compile on mac here
#include <unordered_map>

uint32_t nt4(char);

uint64_t hash64(uint64_t key, const uint64_t& mask);

void test_table();

class KmerHash {
private:
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> lookup;
public:
    static std::pair<uint64_t, uint64_t> kmerhash(const std::string& s, const uint32_t k);
};

#endif
