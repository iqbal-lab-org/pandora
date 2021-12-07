#ifndef PANDORA_DENOVO_RECORD_H
#define PANDORA_DENOVO_RECORD_H

#include <string>
#include <vector>
#include "utils.h"
#include <seqan/align.h>

struct DenovoVariantRecord {
public:
    uint32_t pos;
    std::string ref, alt;
    DenovoVariantRecord(
        uint32_t pos, const std::string& ref, const std::string& alt)
        : pos(pos)
        , ref(ref)
        , alt(alt)
    {
    }
    std::string to_string() const;
    static std::vector<DenovoVariantRecord> get_variants_from_pair_of_sequences(
        const std::string &ref_as_str,
        const std::string &alt_as_str);
};


#endif // PANDORA_DENOVO_RECORD_H
