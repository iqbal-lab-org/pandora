#ifndef PANDORA_SAMPLEINFO_H
#define PANDORA_SAMPLEINFO_H

#include <unordered_map>
#include <string>
#include <vector>
#include <cassert>

template<class STORAGE_TYPE>
class SampleInfo : public std::unordered_map<std::string, std::vector< STORAGE_TYPE > > {
private:
    static const std::vector<std::string> keys_that_can_be_merged;

    void merge_sample_key_from_other_sample_info_into_this (const SampleInfo<STORAGE_TYPE> &other, const std::string &key) {
        if (this->empty() or this->find(key) == this->end() or this->at(key).empty()) {
            return;
        } else if (this->find(key) != this->end() and (other.find(key) == other.end() or other.at(key).empty())) {
            this->erase(key);
        } else if (this->at(key)[0] == other.at(key).at(0)) {
            bool ref = true;
            for (const auto &val : other.at(key)) {
                if (!ref)
                    this->at(key).push_back(val);
                ref = false;
            }
        } else
            this->erase(key);
    }

public:
    void merge_other_sample_info_into_this (const SampleInfo<STORAGE_TYPE> &other) {
        for (const auto &key : keys_that_can_be_merged) {
            merge_sample_key_from_other_sample_info_into_this(other, key);
        }
    }

};

template <class STORAGE_TYPE>
const std::vector<std::string> SampleInfo<STORAGE_TYPE>::keys_that_can_be_merged = {"MEAN_FWD_COVG", "MEAN_REV_COVG",
                                                          "MED_FWD_COVG", "MED_REV_COVG",
                                                          "SUM_FWD_COVG", "SUM_REV_COVG",
                                                          "LIKELIHOOD", "GT_CONF", "GAPS"};



template<class STORAGE_TYPE>
class SamplesInfos : public std::vector<SampleInfo<STORAGE_TYPE>> {
public:
    inline void push_back_several_empty_sample_infos (size_t amount) {
        for (size_t index = 0; index < amount; ++index)
            this->push_back(SampleInfo<STORAGE_TYPE>());
    }

    void merge_other_samples_infos_into_this(const SamplesInfos<STORAGE_TYPE> &other) {
        bool samples_infos_to_be_merged_have_the_same_size = this->size() == other.size();
        assert(samples_infos_to_be_merged_have_the_same_size);

        for (size_t sample_index = 0; sample_index < this->size(); ++sample_index) {
            (*this)[sample_index].merge_other_sample_info_into_this(other[sample_index]);
        }
    }

};

#endif //PANDORA_SAMPLEINFO_H
