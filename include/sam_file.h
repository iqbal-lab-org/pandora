#ifndef PANDORA_SAM_FILE_H
#define PANDORA_SAM_FILE_H

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "forward_declarations.h"
#include "generic_file.h"
#include "seq.h"
#include "localPRG.h"

class Cigar {
private:
    std::vector<char> ops;
    std::vector<uint32_t> lengths;
public:
    Cigar() = default;
    Cigar(const Cigar& other) = default;
    Cigar(Cigar&& other) = default;
    Cigar& operator=(const Cigar& other) = default;
    Cigar& operator=(Cigar&& other) = default;
    virtual ~Cigar() = default;

    inline void add_entry(char op, uint32_t length) {
        ops.emplace_back(op);
        lengths.emplace_back(length);
    }

    inline uint32_t number_of_mismatches() const {
        std::vector<uint32_t> mismatches(lengths.size(), 0);
        std::transform (ops.begin(), ops.end(), lengths.begin(), mismatches.begin(),
            [](char op, uint32_t length) -> uint32_t {
                const bool is_mismatch = op=='X';
                return is_mismatch ? length : 0;
            }
        );
        return std::accumulate(mismatches.begin(), mismatches.end(), (uint32_t)0);
    }

    inline std::string to_string() const {
        std::stringstream ss;
        for(size_t i=0; i<ops.size(); ++i) {
            ss << lengths[i] << ops[i];
        }
        return ss.str();
    }

    inline friend std::ostream& operator<<(std::ostream& out, const Cigar& cigar) {
        out << cigar.to_string();
        return out;
    }
};

class SAMFile : public GenericFile {
private:
    std::vector<bool> get_mapped_positions_bitset(const Seq &seq, const Hits &cluster) const;
    Cigar get_cigar(const std::vector<bool> &mapped_positions_bitset) const;
    std::string get_segment_sequence(const Seq &seq,
                                     const std::vector<bool> &mapped_positions_bitset) const;
    const std::vector<std::shared_ptr<LocalPRG>>& prgs;
public:
    SAMFile(const fs::path &filepath,
            // just to convert prg IDs to prg names
            const std::vector<std::shared_ptr<LocalPRG>>& prgs);
    void write_sam_record_from_hit_cluster(
        const Seq &seq, const MinimizerHitClusters &clusters);

private:
    // some constants
    const static uint32_t flank_size = 50;
};

#endif // PANDORA_SAM_FILE_H
