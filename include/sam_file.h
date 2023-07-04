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

class SAMFile {
private:
    const fs::path filepath;
    const fs::path tmp_filepath;
    GenericFile* tmp_sam_file;
    const std::vector<std::string> &prg_names;
    const std::vector<uint32_t> &prg_lengths;
    std::set<uint32_t> prg_ids_that_we_mapped_to;
    const uint32_t flank_size;

    std::vector<bool> get_mapped_positions_bitset(const Seq &seq, const MinimizerHits &cluster) const;
    Cigar get_cigar(const std::vector<bool> &mapped_positions_bitset) const;
    std::string get_segment_sequence(const Seq &seq,
                                     const std::vector<bool> &mapped_positions_bitset) const;
    std::string get_header() const;
    void create_final_sam_file();
public:
    SAMFile(const fs::path &filepath, const std::vector<std::string> &prg_names,
            const std::vector<uint32_t> &prg_lengths, const uint32_t flank_size)
        :filepath(filepath), tmp_filepath(filepath.string() + ".tmp"),
         tmp_sam_file(new GenericFile(tmp_filepath)), prg_names(prg_names),
         prg_lengths(prg_lengths), flank_size(flank_size){}
    virtual ~SAMFile();

    // Note: these two methods could be one, e.g. write_sam_record() could generate the
    // SAM string and output it. The reason they are separate is solely for performance
    // reasons. Each thread can get the SAM record description with get_sam_record_from_hit_cluster()
    // first, and then synchronise just to write them to disk
    std::string get_sam_record_from_hit_cluster(const Seq &seq, const MinimizerHitClusters &clusters);
    void write_sam_record(const std::string &sam_record_as_str) {
        (*tmp_sam_file) << sam_record_as_str;
    }

};

#endif // PANDORA_SAM_FILE_H
