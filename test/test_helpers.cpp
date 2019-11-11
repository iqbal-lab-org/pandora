#include "test_helpers.h"

GenotypingOptions default_genotyping_options({1,1,1,1,1,1,1,1,1,1}, 0.01, 0, 0, 0, 0, 0, 0, false);

VCF create_VCF_with_default_parameters() {
    VCF vcf(&default_genotyping_options);
    vcf.add_samples({"sample_1"});
    return vcf;
}