#ifndef PANDORA_SAMPLEINFO_H
#define PANDORA_SAMPLEINFO_H

#include <unordered_map>
#include <string>
#include <vector>

template<class STORAGE_TYPE>
class SampleInfo : public std::unordered_map<std::string, std::vector< STORAGE_TYPE > > {
};

template<class STORAGE_TYPE>
class SamplesInfos : public std::vector<SampleInfo<STORAGE_TYPE>> {
public:
    inline void push_back_several_empty_sample_infos (size_t amount) {
        for (size_t index = 0; index < amount; ++index)
            this->push_back(SampleInfo<STORAGE_TYPE>());
    }
};

#endif //PANDORA_SAMPLEINFO_H
