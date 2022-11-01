#ifndef PANDORA_TEST_HELPERS_H
#define PANDORA_TEST_HELPERS_H

#include <gmock/gmock.h>
#include "vcf.h"

using ::testing::HasSubstr;

extern GenotypingOptions default_genotyping_options;

VCF create_VCF_with_default_parameters(size_t nb_of_samples = 1);

// Adapted from https://stackoverflow.com/a/39578934
#define ASSERT_EXCEPTION(TRY_BLOCK, EXCEPTION_TYPE, MESSAGE)                           \
    try {                                                                              \
        {                                                                              \
            TRY_BLOCK;                                                                 \
        }                                                                              \
        FAIL() << "exception '" << MESSAGE << "' not thrown at all!";                  \
    } catch (const EXCEPTION_TYPE& e) {                                                \
        EXPECT_THAT(e.what(), HasSubstr(MESSAGE))                                      \
            << " exception message is incorrect. Expected the following "              \
               "message:\n\n"                                                          \
            << MESSAGE << "\n";                                                        \
    } catch (...) {                                                                    \
        FAIL() << "exception '" << MESSAGE << "' not thrown with expected type '"      \
               << #EXCEPTION_TYPE << "'!";                                             \
    }

#endif // PANDORA_TEST_HELPERS_H
