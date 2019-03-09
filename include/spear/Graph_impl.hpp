#ifndef SPEAR_GRAPH_IMPL_HPP
#define SPEAR_GRAPH_IMPL_HPP

#include "spear/Molecule.hpp"

namespace Spear {

/******************************************************************************
* Implemenation of
* BondEdge
******************************************************************************/

inline AtomVertex BondEdge::source() const {
    return AtomVertex(br_, boost::source(index_, br_->graph()));
}

inline AtomVertex BondEdge::target() const {
    return AtomVertex(br_, boost::target(index_, br_->graph()));
}

inline size_t BondEdge::order() const {
    auto s = source();
    auto t = target();
    return br_->frame().topology().bond_order(s, t);
}

/******************************************************************************
 * Implemenation of
 * AtomVertex::NeighborIterator
 ******************************************************************************/

inline bool AtomVertex::NeighborIterator::operator==(const NeighborIterator& rhs) const {
    return rhs.br_ == this->br_ && rhs.index_ == this->index_;
}

inline bool AtomVertex::NeighborIterator::operator!=(const NeighborIterator& rhs) const {
    return !(*this == rhs);
}

inline AtomVertex::NeighborIterator& AtomVertex::NeighborIterator::operator++() {
    index_++;
}

inline AtomVertex::NeighborIterator AtomVertex::NeighborIterator::operator++(int) {
    NeighborIterator tmp (*this);
    ++(*this);
    return tmp;
}

inline AtomVertex AtomVertex::NeighborIterator::operator* () const {
        return {br_, *index_};
}

/******************************************************************************
* Implemenation of
* AtomVertex::EdgeIterator
******************************************************************************/

inline bool AtomVertex::BondIterator::operator==(const BondIterator& rhs) const {
   return rhs.br_ == this->br_ && rhs.bindex_ == this->bindex_;
}

inline bool AtomVertex::BondIterator::operator!=(const BondIterator& rhs) const {
   return !(*this == rhs);
}

inline AtomVertex::BondIterator& AtomVertex::BondIterator::operator++() {
   bindex_++;
}

inline AtomVertex::BondIterator AtomVertex::BondIterator::operator++(int) {
   BondIterator tmp (*this);
   ++(*this);
   return tmp;
}

inline BondEdge AtomVertex::BondIterator::operator* () const {
       return {br_, *bindex_};
}

/******************************************************************************
 * Implemenation of
 * AtomVertex
 ******************************************************************************/

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
    return static_cast<size_t>(std::distance(tmp.begin(), tmp.end()));
}

inline AtomVertex::Neighbors AtomVertex::neighbors() const {
    return {boost::adjacent_vertices(index_, br_->graph()), br_};
}

inline AtomVertex::Bonds AtomVertex::bonds() const {
    return {boost::out_edges(index_, br_->graph()), br_};
}

inline bool AtomVertex::is_aromatic() const {
    auto types = br_->get_default_atomtype();
    return types->is_aromatic(index_);
}

inline size_t AtomVertex::implicit_hydrogens() const {
    return neighbor_count() < expected_bonds()?
           expected_bonds() - neighbor_count() : 0;
}

inline size_t AtomVertex::explicit_hydrogens() const {
    size_t count = 0;
    auto& mol = *br_;
    for (auto neighbor : this->neighbors()) {
        if (mol[neighbor].atomic_number() == 1) {
            ++count;
        }
    }
    return count;
}

inline size_t AtomVertex::total_hydrogens() const {
    return explicit_hydrogens() + implicit_hydrogens();
}

inline bool AtomVertex::operator==(const AtomVertex& rhs) const {
    return (br_ == rhs.br_ && index_ == rhs.index_);
}

inline AtomVertex AtomVertex::operator[](size_t i) const {
    auto tmp = neighbors().begin();
    std::advance(tmp, i);
    return AtomVertex(br_, *(tmp));
}

inline AtomVertex::operator size_t() const {
    return index_;
}

}

#endif
