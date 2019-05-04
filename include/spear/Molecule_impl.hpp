#ifndef SPEAR_MOLECULE_IMPL_HPP
#define SPEAR_MOLECULE_IMPL_HPP

#include "spear/Molecule.hpp"
#include "spear/Graph_impl.hpp"

namespace Spear {

/******************************************************************************
 * Implemenation of
 * Molecule::iterator
 ******************************************************************************/

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator==(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ == rhs.index_);
}

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator!=(const iterator& rhs) const {
    return !(*this == rhs);
}

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator>=(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ >= rhs.index_);
}

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator<=(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ <= rhs.index_);
}

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator>(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ > rhs.index_);
}

template<typename ret_type, typename sto_type>
inline bool Molecule::iterator<ret_type, sto_type>::operator<(const iterator& rhs) const {
    return (base_mol_ == rhs.base_mol_ && index_ < rhs.index_);
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>&
Molecule::iterator<ret_type, sto_type>::operator++() {
    ++index_;
    return *this;
}   

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>
Molecule::iterator<ret_type, sto_type>::operator++(int) {
    iterator<ret_type, sto_type> tmp (*this);
    ++(*this);
    return tmp;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>&
Molecule::iterator<ret_type, sto_type>::operator--() {
    --index_;
    return *this;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>
Molecule::iterator<ret_type, sto_type>::operator--(int) {
    iterator<ret_type, sto_type> tmp (*this);
    --(*this);
    return tmp;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>&
Molecule::iterator<ret_type, sto_type>::operator+=(difference_type i) {
    std::advance(index_, i);
    return *this;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>&
Molecule::iterator<ret_type, sto_type>::operator-=(difference_type i) {
    index_ -= i;
    return *this;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>
Molecule::iterator<ret_type, sto_type>::operator+(difference_type i) const {
    iterator tmp(*this);
    tmp += i;
    return tmp;
}

template<typename ret_type, typename sto_type>
inline Molecule::iterator<ret_type, sto_type>
Molecule::iterator<ret_type, sto_type>::operator-(difference_type i) const {
    iterator tmp(*this);
    tmp -= i;
    return tmp;
}

template<typename ret_type, typename sto_type>
inline typename Molecule::iterator<ret_type, sto_type>::difference_type
Molecule::iterator<ret_type, sto_type>::operator-(const iterator<ret_type, sto_type>& rhs) const {
    if (base_mol_ != rhs.base_mol_) {
        throw std::logic_error("Cannot take the difference of unrelated iterators");
    }

    return
    index_ >= rhs.index_ ?  static_cast<difference_type>(index_ - rhs.index_)
                         : -static_cast<difference_type>(rhs.index_ - index_);
}

template<typename ret_type, typename sto_type>
inline typename Molecule::iterator<ret_type, sto_type>::value_type
Molecule::iterator<ret_type, sto_type>::operator* () const {
    return ret_type(base_mol_, *index_);
}

template<typename ret_type, typename sto_type>
inline typename Molecule::iterator<ret_type, sto_type>::value_type
Molecule::iterator<ret_type, sto_type>::operator[](difference_type rhs) const {
    iterator tmp(*this);
    tmp += rhs;
    return *tmp;
}

/******************************************************************************
 * Implemenation of
 * Molecule
 ******************************************************************************/

inline BondEdge Molecule::bond(size_t idx1, size_t idx2) {
    auto edge = boost::edge(idx1, idx2, graph_);
    if (!edge.second) {
        throw std::runtime_error("No bond between: " + std::to_string(idx1) +
                                 " and " + std::to_string(idx2));
    }
    return {this, edge.first};
}

template<typename atomtype, typename... args>
inline std::string Molecule::add_atomtype(args... additional) {
    auto typed_atoms = new atomtype(*this, additional...);
    auto name = typed_atoms->name();
    atom_types_[name] = std::unique_ptr<atomtype>(typed_atoms);

    if (default_atomtype_.empty()) {
        default_atomtype_ = name;
    }

    return name;
}

inline const AtomType* Molecule::atomtype(const std::string& name) const {
    auto types = name.empty() ? atom_types_.find(default_atomtype_)
                              : atom_types_.find(name);
    if (types == atom_types_.end()) {
        return nullptr;
    }
    return types->second.get();
}

inline void Molecule::set_default_atomtype(const std::string& name) {
    auto types = atom_types_.find(name);
    if (types == atom_types_.end()) {
        throw std::runtime_error("Atom type " + name + " is not availible.");
    }
    default_atomtype_ = name;
}

inline size_t Molecule::size() const {
    return topology_.size();
}

inline AtomVertex Molecule::operator[](size_t index) const {
    if (index >= size()) {
        throw std::out_of_range("Index given to Molecule::operator[] is too large.");
    }
    return AtomVertex(this, index);
}

inline Molecule::iterator<AtomVertex, VertexIterator> Molecule::begin() const {
    auto vert_iters = boost::vertices(graph_);
    return iterator<AtomVertex, VertexIterator>(this, vert_iters.first);
}

inline Molecule::iterator<AtomVertex, VertexIterator> Molecule::end() const {
    auto vert_iters = boost::vertices(graph_);
    return iterator<AtomVertex, VertexIterator>(this, vert_iters.second);
}

inline Molecule::iterator<AtomVertex, VertexIterator> Molecule::cbegin() const {
    auto vert_iters = boost::vertices(graph_);
    return iterator<AtomVertex, VertexIterator>(this, vert_iters.first);
}

inline Molecule::iterator<AtomVertex, VertexIterator> Molecule::cend() const {
    auto vert_iters = boost::vertices(graph_);
    return iterator<AtomVertex, VertexIterator>(this, vert_iters.second);
}

inline Molecule::AllBonds Molecule::bonds() const {
    return AllBonds(this, boost::edges(graph_));
}

}

#endif
