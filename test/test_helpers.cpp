#include "test_helpers.h"

GenotypingOptions default_genotyping_options(
    { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 0.01, 0, 0, 0, 0, 0, 0, false);

VCF create_VCF_with_default_parameters(size_t nb_of_samples)
{
    VCF vcf(&default_genotyping_options);

    std::vector<std::string> sample_names;
    for (size_t sample_index = 1; sample_index <= nb_of_samples; ++sample_index) {
        sample_names.push_back(std::string("sample_") + std::to_string(sample_index));
    }
    vcf.add_samples(sample_names);

    return vcf;
}