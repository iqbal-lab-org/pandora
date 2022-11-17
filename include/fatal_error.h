#ifndef PANDORA_FATAL_ERROR_H
#define PANDORA_FATAL_ERROR_H

#include <sstream>

template <typename T> std::string to_string(const T& element)
{
    std::stringstream ss;
    ss << element;
    return ss.str();
}

// From https://stackoverflow.com/a/21806609
template <typename... Args> std::string stringer(Args const&... args)
{
    std::string result;
    using ::to_string;
    using std::to_string;
    int unpack[] { 0, (result += to_string(args), 0)... };
    static_cast<void>(unpack);
    return result;
}

class FatalRuntimeError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

template <typename... Args> void fatal_error(Args const&... args)
{
    throw FatalRuntimeError(stringer("[FATAL ERROR]: ", args..., "\nAborting..."));
}

#endif // PANDORA_FATAL_ERROR_H
