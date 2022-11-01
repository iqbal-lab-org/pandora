#include "denovo_discovery/denovo_record.h"

std::string DenovoVariantRecord::to_string() const
{
    std::stringstream ss;
    ss << pos << "\t" << ref << "\t" << alt;
    return ss.str();
}

void add_and_reset_variant(std::vector<DenovoVariantRecord> &denovo_variants,
    int64_t &variant_pos, std::string &variant_ref, std::string &variant_alt) {
    const bool should_add_variant = variant_pos != -1;
    if (should_add_variant) {
        denovo_variants.emplace_back((uint32_t)variant_pos, variant_ref,
            variant_alt);
        variant_pos = -1;
        variant_ref = variant_alt = "";
    }
}

std::vector<DenovoVariantRecord> DenovoVariantRecord::get_variants_from_pair_of_sequences(
    const std::string &ref_as_str,
    const std::string &alt_as_str)
{
    using namespace seqan;
    typedef String<char> TSequence; // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign; // align type
    typedef Row<TAlign>::Type TRow; // gapped sequence type

    TSequence ref(ref_as_str);
    TSequence alt(alt_as_str);

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), ref);
    assignSource(row(align, 1), alt);

    TRow& ref_row = row(align, 0);
    TRow& alt_row = row(align, 1);

    globalAlignment(align, Score<int, Simple>(2, -1, -2, -4), AffineGaps());
    bool start_of_variant = true;
    std::vector<DenovoVariantRecord> denovo_variants;
    std::string variant_ref, variant_alt;
    int64_t variant_pos=-1;
    for (size_t alignment_index = 0; alignment_index < length(ref_row);
         ++alignment_index) {
        char ref_base = (char)(ref_row[alignment_index]);
        char alt_base = (char)(alt_row[alignment_index]);
        const bool is_variant_position = ref_base != alt_base;
        if (is_variant_position) {
            if (start_of_variant) {
                variant_pos = toSourcePosition(ref_row, alignment_index) + 1;
            }
            if (ref_base != '-') variant_ref += ref_base;
            if (alt_base != '-') variant_alt += alt_base;
            start_of_variant = false;
        } else {
            add_and_reset_variant(denovo_variants, variant_pos, variant_ref, variant_alt);
            start_of_variant = true;
        }
    }
    add_and_reset_variant(denovo_variants, variant_pos, variant_ref, variant_alt);

    return denovo_variants;
}
