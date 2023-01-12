#ifndef PANDORA_CLI_HELPERS_H
#define PANDORA_CLI_HELPERS_H

#include "CLI11.hpp"

class PandoraIndexValidator : public CLI::Validator {
public:
    PandoraIndexValidator() : Validator("") {
        func_ = [](std::string &filename) {
            const bool index_has_the_correct_extension =
                filename.rfind(".panidx.zip") == filename.size()-11;
            if (index_has_the_correct_extension) {
                return std::string();
            }
            else {
                return std::string("Index file has the wrong extension, value: ") +
                       filename +
                       " (expected *.panidx.zip)";
            }
        };
    }
};


#endif // PANDORA_CLI_HELPERS_H
