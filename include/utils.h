#ifndef __UTILS_H_INCLUDED__   // if utils.h hasn't been included yet...
#define __UTILS_H_INCLUDED__

template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};
#endif
