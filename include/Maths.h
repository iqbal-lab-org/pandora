#ifndef PANDORA_MATHS_H
#define PANDORA_MATHS_H

#include <iterator>
#include <functional>
#include <numeric>
#include <algorithm>
#include <vector>
#include "gtest/gtest.h"
#include <cmath>

class Maths {
public:
    template<class Iterator>
    inline static typename std::iterator_traits<Iterator>::value_type sum(Iterator begin, Iterator end) {
        typedef typename std::iterator_traits<Iterator>::value_type value_type;

        value_type default_value = value_type();
        return std::accumulate(begin, end, default_value);
    }

    template<class Iterator>
    inline static typename std::iterator_traits<Iterator>::value_type mean(Iterator begin, Iterator end) {
        typedef typename std::iterator_traits<Iterator>::difference_type difference_type;

        difference_type number_of_elements = std::distance(begin, end);
        bool no_elements_in_container = number_of_elements == 0;
        if (no_elements_in_container) {
            throw std::runtime_error("no_elements_in_container at Maths::mean()");
        }

        return Maths::sum(begin, end) / number_of_elements;
    }


    template<class Iterator>
    inline static typename std::iterator_traits<Iterator>::value_type median(Iterator begin, Iterator end) {
        // TODO: performance improvement: copy pointers to the vector, not the elements themselves
        typedef typename std::iterator_traits<Iterator>::value_type value_type;

        std::vector<value_type> copy_vector(begin, end);
        if (copy_vector.empty()) {
            throw std::runtime_error("no_elements_in_container at Maths::median()");
        }

        std::sort(copy_vector.begin(), copy_vector.end());
        if (copy_vector.size() % 2 == 1) {
            size_t n = (copy_vector.size() + 1) / 2;
            return copy_vector[n - 1];
        } else {
            size_t n1 = (copy_vector.size() + 1) / 2;
            size_t n2 = (copy_vector.size() - 1) / 2;
            return (copy_vector[n1] + copy_vector[n2]) / 2;
        }
    }


    template<class Iterator>
    inline static typename std::iterator_traits<Iterator>::value_type mode(Iterator begin, Iterator end) {
        // TODO: performance improvement: copy pointers to the vector, not the elements themselves
        typedef typename std::iterator_traits<Iterator>::value_type value_type;

        std::vector<value_type> copy_vector(begin, end);
        if (copy_vector.empty()) {
            throw std::runtime_error("no_elements_in_container at Maths::mode()");
        }

        std::sort(copy_vector.begin(), copy_vector.end());
        uint32_t counter = 1;
        uint32_t max_count = 1;
        value_type most_common_value = copy_vector[0];
        value_type last_value = copy_vector[0];
        for_each(copy_vector.begin()+1, copy_vector.end(), [&](const value_type &value) {
            if (value == last_value)
                counter++;
            else {
                if (counter > max_count) {
                    max_count = counter;
                    most_common_value = last_value;
                }
                counter = 1;
            }
            last_value = value;
        });

        if (counter > max_count) {
            max_count = counter;
            most_common_value = last_value;
        }

        return most_common_value;
    }

    // floating point comparison is not well defined, these functions compare floating pointers based on how google test does it
    // based on https://stackoverflow.com/a/3423299/5264075
    template <class FloatOrDouble>
    inline static bool equals(FloatOrDouble left, FloatOrDouble right) {
        const testing::internal::FloatingPoint<FloatOrDouble> lhs(left), rhs(right);
        return lhs.AlmostEquals(rhs);
    }

    inline static double logfactorial (uint32_t n) {
        double logfactorial = 0;
        for (uint32_t i = 1; i <= n; ++i) {
            logfactorial += std::log(i);
        }
        return logfactorial;
    }

    template<class Iterator>
    inline static typename std::iterator_traits<Iterator>::difference_type arg_max (Iterator begin, Iterator end) {
        auto max_element_iterator = std::max_element(begin, end);
        auto index_of_max_element = std::distance(begin, max_element_iterator);
        return index_of_max_element;
    }

};


#endif //PANDORA_MATHS_H
