#include "denovo_discovery/denovo_record.h"

std::string DenovoVariantRecord::to_string() const
{
    std::string ref_to_print = remove_spaces_from_string(ref);
    std::string alt_to_print = remove_spaces_from_string(alt);
    std::stringstream ss;
    ss << pos << "\t" << ref_to_print << "\t" << alt_to_print;
    return ss.str();
}

std::vector<DenovoVariantRecord> DenovoVariantRecord::get_variants_from_pair_of_sequences(
    const std::string &ref_as_str,
    const std::string &alt_as_str)
{
    using namespace seqan;
    typedef String<char> TSequence; // sequence type
    typedef Align<TSequence, ArrayGaps> TAlign; // align type
    typedef Row<TAlign>::Type TRow; // gapped sequence type

    // TODO: this can be further optimised by aligning the candidate region sequence
    // with the candidate region alt (not the whole ML sequence with the whole ML alt)
    TSequence ref(ref_as_str);
    TSequence alt(alt_as_str);

    TAlign align;
    resize(rows(align), 2);
    assignSource(row(align, 0), ref);
    assignSource(row(align, 1), alt);

    TRow& ref_row = row(align, 0);
    TRow& alt_row = row(align, 1);

    globalAlignment(align, Score<int, Simple>(2, -1, -2, -4), AffineGaps());
    bool append_to_previous = false;
    std::vector<DenovoVariantRecord> denovo_variants;
    for (size_t alignment_index = 0; alignment_index < length(ref_row);
         ++alignment_index) {
        char ref_base = (char)(ref_row[alignment_index]);
        char alt_base = (char)(alt_row[alignment_index]);
        const bool is_variant_position = ref_base != alt_base;
        if (is_variant_position) {
            if (append_to_previous) {
                denovo_variants.back().ref += ref_base;
                denovo_variants.back().alt += alt_base;
            } else {
                DenovoVariantRecord denovo_variant(
                    toSourcePosition(ref_row, alignment_index) + 1,
                    std::string(1, ref_base), std::string(1, alt_base));
                denovo_variants.push_back(denovo_variant);
            }
            append_to_previous = true;
        } else {
            append_to_previous = false;
        }
    }

    return denovo_variants;
}
