#ifndef SPEAR_MOLECULE_IMPL_HPP
#define SPEAR_MOLECULE_IMPL_HPP

#include "spear/Molecule.hpp"
#include "spear/Graph_impl.hpp"

namespace Spear {

/******************************************************************************
 * Implemenation of
 * Molecule::iterator
 ******************************************************************************/

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

/******************************************************************************
 * Implemenation of
 * Molecule
 ******************************************************************************/

template<class atomtype, typename... args>
inline std::string Molecule::add_atomtype(args... additional) {
    auto typed_atoms = new atomtype(*this, additional...);
    auto name = typed_atoms->name();
    atom_types_[name] = std::unique_ptr<atomtype>(typed_atoms);

    if (default_atomtype_.empty()) {
        default_atomtype_ = name;
    }

    return name;
}

inline optional<const AtomType*> Molecule::get_atomtype(const std::string& name) const {
    auto types = atom_types_.find(name);
    if (types == atom_types_.end()) {
        return chemfiles::nullopt;
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

inline const AtomType* Molecule::get_default_atomtype() const {
    auto types = atom_types_.find(default_atomtype_);
    return types->second.get();
}

inline size_t Molecule::size() const {
    return frame_.size();
}

inline AtomVertex Molecule::operator[](size_t index) const {
    if (index >= size()) {
        throw std::out_of_range("Index given to Molecule::operator[] is too large.");
    }
    return AtomVertex(this, index);
}

inline Molecule::iterator Molecule::begin() const {
    return iterator(this, 0);
}

inline Molecule::iterator Molecule::end() const {
    return iterator(this, size());
}

inline Molecule::iterator Molecule::cbegin() const {
    return iterator(this, 0);
}

inline Molecule::iterator Molecule::cend() const {
    return iterator(this, size());
}

}

#endif
