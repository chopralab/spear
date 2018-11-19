// License:

#ifndef SPEAR_MOLECULE_HPP
#define SPEAR_MOLECULE_HPP

#include "spear/exports.hpp"
#include "chemfiles/Frame.hpp"

#include <vector>
#include <boost/graph/undirected_graph.hpp>

namespace Spear {

// Note: If more atom/bond properties are need other than number and order,
// change the lines below:
using VertexProperty = boost::property<boost::vertex_name_t, uint64_t>;
using EdgeProperty = boost::property<boost::edge_name_t, uint64_t>;

using Graph = boost::adjacency_list<boost::setS, boost::vecS,boost::undirectedS,
                                    VertexProperty, EdgeProperty>;

using Traits = boost::graph_traits<Graph>;
using VertexDescriptor = Traits::vertex_descriptor;
using EdgeDescriptor = Traits::edge_descriptor;
using AdjacencyIterator = Traits::adjacency_iterator;
using EdgeIterator = Traits::edge_iterator;
using AdjacencyIteratorPair = std::pair<AdjacencyIterator, AdjacencyIterator>;

class Molecule;

class AtomVertex {
public:
    AtomVertex(const Molecule* br, size_t index);

    const std::string& name() const;

    const std::string& type() const;
    
    const chemfiles::Vector3D& position() const;

    size_t atomic_number() const;

    size_t neighbor_count() const;
    
    AdjacencyIteratorPair neighbors() const;

    bool operator==(const AtomVertex& rhs) const;

    AtomVertex operator[](size_t) const;

    operator size_t() const;
    
    AdjacencyIterator begin() const;
    
    AdjacencyIterator end() const;

private:
    size_t index_;
    const Molecule* br_;
};

class SPEAR_EXPORT Molecule {
public:

    class iterator : public std::iterator<std::random_access_iterator_tag, AtomVertex> {
    public:

        using difference_type =
            std::iterator<std::random_access_iterator_tag, AtomVertex>::difference_type;

        iterator(const Molecule* base_mol, size_t i = 0) : index_(i), base_mol_(base_mol) {
        }

        bool operator==(const iterator& rhs) const;

        bool operator!=(const iterator& rhs) const;

        bool operator>=(const iterator& rhs) const;

        bool operator<=(const iterator& rhs) const;

        bool operator>(const iterator& rhs) const;

        bool operator<(const iterator& rhs) const;

        iterator& operator++();

        iterator operator++(int);

        iterator& operator--();

        iterator operator--(int);

        iterator& operator+=(difference_type i);

        iterator& operator-=(difference_type i);

        iterator operator+(difference_type i) const;

        iterator operator-(difference_type i) const;

        difference_type operator-(const iterator& rhs) const;

        AtomVertex operator* () const;

        AtomVertex operator[](difference_type rhs) const;

    private:
        size_t index_;
        const Molecule* base_mol_;
    };

    explicit Molecule(const chemfiles::Frame& frame);

    const Graph& graph() const {
        return graph_;
    }

    const chemfiles::Frame& frame() const {
        return frame_;
    }

    const std::set<std::set<size_t>> rings() const;

    std::vector<EdgeDescriptor> get_bonds_in(const std::set<size_t>& atoms) const;

    size_t size() const {
        return frame_.size();
    }

    AtomVertex operator[](size_t index) const {
        if (index >= size()) {
            throw std::out_of_range("Index given to Molecule::operator[] is too large.");
        }
        return AtomVertex(this, index);
    }

    iterator begin() const {
        return iterator(this, 0);
    }

    iterator end() const {
        return iterator(this, size());
    }

    iterator cbegin() const {
        return iterator(this, 0);
    }

    iterator cend() const {
        return iterator(this, size());
    }

private:
    chemfiles::Frame frame_;
    
    Graph graph_;

    std::vector<size_t> atom_types_;
};

inline AtomVertex::AtomVertex(const Molecule* br, size_t index) :
    index_(index), br_(br) {
}

inline const std::string& AtomVertex::name() const {
    return br_->frame()[index_].name();
}

inline const std::string& AtomVertex::type() const {
    return br_->frame()[index_].type();
}

inline const chemfiles::Vector3D& AtomVertex::position() const {
    return br_->frame().positions()[index_];
}

inline uint64_t AtomVertex::atomic_number() const {
    return *(br_->frame()[index_].atomic_number());
}

inline size_t AtomVertex::neighbor_count() const {
    auto tmp = neighbors();
    return static_cast<size_t>(std::distance(tmp.first, tmp.second));
}

inline AdjacencyIteratorPair AtomVertex::neighbors() const {
    return boost::adjacent_vertices(index_, br_->graph());
}

inline bool AtomVertex::operator==(const AtomVertex& rhs) const {
    return (br_ == rhs.br_ && index_ == rhs.index_);
}

inline AtomVertex AtomVertex::operator[](size_t i) const {
    auto tmp = begin();
    std::advance(tmp, i);
    return AtomVertex(br_, *(tmp));
}

inline AtomVertex::operator size_t() const {
    return index_;
}

inline AdjacencyIterator AtomVertex::begin() const {
    return boost::adjacent_vertices(index_, br_->graph()).first;
}

inline AdjacencyIterator AtomVertex::end() const {
    return boost::adjacent_vertices(index_, br_->graph()).second;
}

inline bool Molecule::iterator::operator==(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ == rhs.index_);
}

inline bool Molecule::iterator::operator!=(const iterator& rhs) const {
    return !(*this == rhs);
}

inline bool Molecule::iterator::operator>=(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ >= rhs.index_);
}

inline bool Molecule::iterator::operator<=(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ <= rhs.index_);;
}

inline bool Molecule::iterator::operator>(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ > rhs.index_);
}

inline bool Molecule::iterator::operator<(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ < rhs.index_);;
}

inline Molecule::iterator& Molecule::iterator::operator++() {
    if(index_ != base_mol_->size())
        ++index_;
    return *this;
}   

inline Molecule::iterator Molecule::iterator::operator++(int) {
    iterator tmp (*this);
    ++(*this);
    return tmp;
}

inline Molecule::iterator& Molecule::iterator::operator--() {
    if(index_ != 0)
        --index_;
    return *this;
}

inline Molecule::iterator Molecule::iterator::operator--(int) {
    iterator tmp (*this);
    --(*this);
    return tmp;
}

inline Molecule::iterator& Molecule::iterator::operator+=(difference_type i) {
    if (i < 0) {
        return (*this) -= (-i);
    }

    auto si = static_cast<size_t>(i);
    index_ = index_ + si > base_mol_->size() ? base_mol_->size() : index_ + si;
    return *this;
}

inline Molecule::iterator& Molecule::iterator::operator-=(difference_type i) {
    if (i < 0) {
        return (*this) += (-i);
    }

    auto si = static_cast<size_t>(i);
    index_ = si > index_ ? 0 : index_ - si;
    return *this;
}

inline Molecule::iterator Molecule::iterator::operator+(difference_type i) const {
    iterator tmp(*this);
    tmp += i;
    return tmp;
}

inline Molecule::iterator Molecule::iterator::operator-(difference_type i) const {
    iterator tmp(*this);
    tmp -= i;
    return tmp;
}

inline Molecule::iterator::difference_type Molecule::iterator::operator-(const iterator& rhs) const {
    if (base_mol_ != rhs.base_mol_) {
        throw std::logic_error("Cannot take the difference of unrelated iterators");
    }

    return
    index_ >= rhs.index_ ?  static_cast<difference_type>(index_ - rhs.index_)
                         : -static_cast<difference_type>(rhs.index_ - index_);
                              
}

inline AtomVertex Molecule::iterator::operator* () const {
    if(index_ == base_mol_->size()) {
        throw std::range_error("Cannot dereference an end iterator!");
    }
    return AtomVertex(base_mol_, index_);
}

inline AtomVertex Molecule::iterator::operator[](difference_type rhs) const {
    iterator tmp(*this);
    tmp += rhs;
    return *tmp;
}

}

#endif
