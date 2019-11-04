#include "test_helpers.h"

GenotypingOptions default_genotyping_options({10,10,10,10,10,10,10,10,10,10}, 0.01, 0, 0, 0, 0, 0, 0, false);

VCF create_VCF_with_default_parameters() {
    return VCF(&default_genotyping_options);
}