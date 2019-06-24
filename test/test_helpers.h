//
// Created by leandro on 14/06/19.
//

#ifndef PANDORA_TEST_HELPERS_H
#define PANDORA_TEST_HELPERS_H

#include <algorithm>

template <typename T,
          template <typename ELEM_TYPE, typename = std::allocator<ELEM_TYPE> > class CONT_TYPE>
bool equal_containers(const CONT_TYPE<T> &lhs, const CONT_TYPE<T> &rhs) {
    return lhs.size()==rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

template <typename T,
          template <typename ELEM_TYPE, typename = std::allocator<ELEM_TYPE> > class CONT_TYPE,
          class BinaryPredicate>
bool equal_containers(const CONT_TYPE<T> &lhs, const CONT_TYPE<T> &rhs, const BinaryPredicate &pred) {
    return lhs.size()==rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin(), pred);
}


#endif //PANDORA_TEST_HELPERS_H
