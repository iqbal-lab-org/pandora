#ifndef __PATH_H_INCLUDED__   // if path.h hasn't been included yet...
#define __PATH_H_INCLUDED__

#include <deque>
#include <vector>
#include <iostream>
#include <cstdint>
#include "interval.h"
#include "prg/ns.cpp"
#include <memory>

class LocalPRG;
class LocalNode;
typedef std::shared_ptr<LocalNode> LocalNodePtr; //TODO: this is replicated from localnode.h and I really don't like it, fix

class prg::Path {
private:
    std::vector<Interval> path; //the interval path - we control acess to this variable now


    //variables for memoization:
    bool isMemoized; //flag that indicated if the first memoization was already done
    std::vector<LocalNodePtr> memoizedLocalNodePath; //the memoized local node path
    bool memoizedDirty; //was this path modified and thus memoizedLocalNodePath needs to be recomputed?
    uint32_t localPRGIdOfMemoizedLocalNodePath; //just to make sure we don't memoize one interval for one localPRG and return it to a different localPRg

public:
    //constructors
    Path()
            : path{}, isMemoized{false}, memoizedLocalNodePath{}, memoizedDirty{true},
              localPRGIdOfMemoizedLocalNodePath{0} {} //default constructor
    Path(const Path &other) = default; //copy default constructor
    Path(Path &&other) { //move constructor
        *this = std::move(other);
    }

    //assignment operators
    Path &operator=(const Path &other) = default; //copy assignment operator
    Path &operator=(Path &&other) {
        if (this != &other) {
            this->path = std::move(other.path);
            this->isMemoized = std::move(other.isMemoized);
            this->memoizedLocalNodePath = std::move(other.memoizedLocalNodePath);
            this->memoizedDirty = std::move(other.memoizedDirty);
            this->localPRGIdOfMemoizedLocalNodePath = std::move(other.localPRGIdOfMemoizedLocalNodePath);
        }
    }

    //destructor
    virtual ~Path() = default;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //DIRTY METHODS - THAT CAN MODIFY path - memoization must be redone
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //delegated std::vector methods
    void push_back(const Interval &interval) {
        memoizedDirty = true;
        path.push_back(interval);
    }

    void clear() {
        memoizedDirty = true;
        path.clear();
    }

    //add all intervals in the given range to the end of the path
    template<class Iterator>
    void insert_to_the_end(const Iterator &begin, const Iterator &end) {
        memoizedDirty = true;
        path.insert(path.end(), begin, end);
    }

    Interval getAndRemoveLastInterval() {
        memoizedDirty = true;
        Interval interval = path.back();
        path.pop_back();
        return interval;
    }


    /////////////////////////////
    //initializers
    //this was not before, but it is the most generic way of intializing (just give any two iterators)
    template<class Iterator>
    void initialize(const Iterator &begin, const Iterator &end, uint32_t reservedSize = 16) {
        memoizedDirty = true;
        path.clear();
        path.reserve(reservedSize);
        insert_to_the_end(begin, end);
    }

    //some other convenience initializers, also because these are spread out in the code
    //initialize with a single interval
    void initialize(const Interval &i) {
        memoizedDirty = true;
        std::vector<Interval> vectorOfInterval = {i};
        initialize(vectorOfInterval.begin(), vectorOfInterval.end());
    }

    //initializes this path with the intervals in the container
    template<class ContainerType>
    void initialize(const ContainerType &container) {
        /* TODO: keep this?
        if (container.empty())
            return;
        */
        memoizedDirty = true;
        initialize(container.begin(), container.end(), container.size() + 5);
    }

    //initializes from a previous path
    void initialize(const prg::Path &p) {
        memoizedDirty = true;
        initialize(p.path);
    }

    //add an interval to the end
    void add_end_interval(const Interval &);
    /////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //DIRTY METHODS - THAT CAN MODIFY path
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CONST METHODS - THAT CANNOT MODIFY path - memoization does not need to be redone in these cases
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //delegated std::vector methods
    bool empty() const { return path.empty(); }

    std::vector<Interval>::const_iterator
    begin() const noexcept { return path.cbegin(); } //iteration on path always use constant iterator for now (no changes to intervals using iterators)
    std::vector<Interval>::const_iterator
    end() const noexcept { return path.cend(); } //iteration on path always use constant iterator for now (no changes to intervals using iterators)
    const Interval &operator[](size_t n) const { return path[n]; }

    const Interval &back() const { return path.back(); }

    size_t size() const { return path.size(); }

    const std::vector<Interval> &getPath() const { return path; }


    //some getters
    uint32_t get_start() const;

    uint32_t get_end() const;

    uint32_t length() const;

    Path subpath(const uint32_t, const uint32_t) const;

    bool is_branching(const Path &) const;

    bool is_subpath(const Path &) const;

    //comparators
    bool operator<(const Path &y) const;

    bool operator==(const Path &y) const;

    bool operator!=(const Path &y) const;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //CONST METHODS - THAT CANNOT MODIFY path
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //MEMOIZATION METHOD
    //the method responsible for memoization itself
    std::vector<LocalNodePtr> nodes_along_path(const LocalPRG &localPrg);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //friend methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend std::ostream &operator<<(std::ostream &out, const Path &p);

    friend std::istream &operator>>(std::istream &in, Path &p);

    friend Path get_union(const Path &, const Path &);

    friend bool equal_except_null_nodes(const Path &, const Path &);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //friend methods
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
};

bool equal_except_null_nodes(const prg::Path &, const prg::Path &);




#endif
