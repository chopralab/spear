#ifndef SPEAR_GRAPH_IMPL_HPP
#define SPEAR_GRAPH_IMPL_HPP

#include "spear/Graph.hpp"
#include "spear/Molecule.hpp"

#include <boost/iterator/distance.hpp>

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

inline AtomVertex BondEdge::other_atom(const AtomVertex& av) const {
    if (source() == av) {
        return source();
    }

    if (target() == av) {
        return target();
    }

    //TODO: throw
    return source();
}

inline std::pair<AtomVertex, AtomVertex> BondEdge::as_pair() const {
    return std::minmax(source(), target());
}

inline size_t BondEdge::index() const {
    auto all_bonds = br_->bonds();

    size_t count = 0;

    for (auto bond : all_bonds) {
        if (as_pair() == bond.as_pair()) {
            break;
        }
        ++count;
    }

    return count;
}

inline Bond::Order BondEdge::order() const {
    return boost::get(boost::edge_name, br_->graph(), index_);
}

/******************************************************************************
 * Implemenation of
 * Neighbors::NeighborIterator
 ******************************************************************************/

inline bool Neighbors::NeighborIterator::operator==(const NeighborIterator& rhs) const {
    return rhs.br_ == this->br_ && rhs.index_ == this->index_;
}

inline bool Neighbors::NeighborIterator::operator!=(const NeighborIterator& rhs) const {
    return !(*this == rhs);
}

inline Neighbors::NeighborIterator& Neighbors::NeighborIterator::operator++() {
    index_++;
	return *this;
}

inline Neighbors::NeighborIterator Neighbors::NeighborIterator::operator++(int) {
    NeighborIterator tmp (*this);
    ++(*this);
    return tmp;
}

inline AtomVertex Neighbors::NeighborIterator::operator* () const {
    return {br_, *index_};
}

/******************************************************************************
* Implemenation of
* Bonds::EdgeIterator
******************************************************************************/

inline bool Bonds::BondIterator::operator==(const BondIterator& rhs) const {
   return rhs.br_ == this->br_ && rhs.bindex_ == this->bindex_;
}

inline bool Bonds::BondIterator::operator!=(const BondIterator& rhs) const {
   return !(*this == rhs);
}

inline Bonds::BondIterator& Bonds::BondIterator::operator++() {
   bindex_++;
   return *this;
}

inline Bonds::BondIterator Bonds::BondIterator::operator++(int) {
   BondIterator tmp (*this);
   ++(*this);
   return tmp;
}

inline BondEdge Bonds::BondIterator::operator* () const {
       return {br_, *bindex_};
}

/******************************************************************************
 * Implemenation of
 * AtomVertex
 ******************************************************************************/

inline AtomVertex::AtomVertex(const Molecule* br, VertexDescriptor index) :
    index_(index), br_(br) {
}

inline const std::string& AtomVertex::name() const {
    return br_->topology()[index_].name();
}

inline const Eigen::Vector3d& AtomVertex::position() const {
    return br_->positions()[index_];
}

inline Element::Symbol AtomVertex::atomic_number() const {
    return boost::get(boost::vertex_name, br_->graph(), index_);
}

inline size_t AtomVertex::degree() const {
    return boost::degree(index_, br_->graph());
}

inline Neighbors AtomVertex::neighbors() const {
    return {boost::adjacent_vertices(index_, br_->graph()), br_};
}

inline Bonds AtomVertex::bonds() const {
    return {boost::out_edges(index_, br_->graph()), br_};
}

inline bool AtomVertex::is_aromatic() const {
    return br_->atomtype()->is_aromatic(index_);
}

inline bool AtomVertex::is_planar() const {
    return br_->atomtype()->is_planar(index_);
}

inline bool AtomVertex::is_non_metal() const {
    auto element = atomic_number();
    if (element <= Element::He ||
       (element >= Element::B  && element <= Element::Ne) ||
       (element >= Element::Si && element <= Element::Ar) ||
       (element >= Element::Ge && element <= Element::Kr) ||
       (element >= Element::Sb && element <= Element::Xe) ){
        return true;
    }

    return false;
}

inline size_t AtomVertex::implicit_hydrogens() const {
    return degree() < expected_bonds()?
           expected_bonds() - degree() : 0;
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

inline Atom::Chirality AtomVertex::chirality() const {
    return Atom::UNSPECIFIED;
}

inline AtomRingMapIteratorPair AtomVertex::rings() const {
    return br_->atom_to_ring_.equal_range(index_);
}

inline AtomRingMapIteratorPair AtomVertex::sssrs() const {
    return br_->atom_to_sssr_.equal_range(index_);
}

inline Hybridization AtomVertex::hybridization() const {
    return br().atomtype()->hybridization(index_);
}

inline double AtomVertex::partial_charge() const {
    return (*br().partial_charge())[index_];
}

inline int8_t AtomVertex::formal_charge() const {
    //TODO: Implement
    return 0;
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
