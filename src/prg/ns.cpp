#include <iostream>


namespace prg {
    class Path;

    Path get_union(const Path &x, const Path &y);

    std::ostream &operator<<(std::ostream &out, const Path &p);

    std::istream &operator>>(std::istream &in, Path &p);
}
