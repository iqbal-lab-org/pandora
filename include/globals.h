#ifndef PANDORA_GLOBALS_H
#define PANDORA_GLOBALS_H

#include <string>

class PandoraGlobals{
public:
    static std::string command_line;
};

#define INDEXING_UPPER_BOUND_DEFAULT 10000000
#define ESTIMATED_INDEX_SIZE_DEFAULT 100000

#endif // PANDORA_GLOBALS_H
