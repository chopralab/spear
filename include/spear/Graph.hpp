#ifndef SPEAR_GRAPH_HPP
#define SPEAR_GRAPH_HPP

#include "spear/Constants.hpp"
#include "spear/Typedefs.hpp"
#include "spear/Rings.hpp"
#include "spear/Geometry.hpp"

namespace Spear {

class Molecule;
class AtomVertex;
class BondEdge;

struct SPEAR_EXPORT Neighbors {

    class SPEAR_EXPORT NeighborIterator {
    public:
        using difference_type = std::ptrdiff_t;

        using iterator_category = std::forward_iterator_tag;
        using value_type = AtomVertex;
        using pointer = AtomVertex*;
        using reference = AtomVertex&;

        NeighborIterator(AdjacencyIterator be, const Molecule* br) :
            index_(be), br_(br) {}

        bool operator==(const NeighborIterator& rhs) const;

        bool operator!=(const NeighborIterator& rhs) const;

        NeighborIterator& operator++();

        NeighborIterator operator++(int);

        AtomVertex operator* () const;
    private:
        AdjacencyIterator index_;
        const Molecule* br_;
    };

    Neighbors(AdjacencyIteratorPair be, const Molecule* br) :
        begin_end_(be), br_(br) {}

    NeighborIterator begin() const {
        return {begin_end_.first, br_};
    }

    NeighborIterator end() const {
        return {begin_end_.second, br_};
    }
private:
    const AdjacencyIteratorPair begin_end_;
    const Molecule* br_;
};

struct SPEAR_EXPORT Bonds {

    class SPEAR_EXPORT BondIterator {
    public:
        using difference_type = std::ptrdiff_t;

        using iterator_category = std::forward_iterator_tag;
        using value_type = AtomVertex;
        using pointer = AtomVertex*;
        using reference = AtomVertex&;

        BondIterator() : br_(nullptr), bindex_() {}

        BondIterator(OutEdgeIterator be, const Molecule* br) :
            br_(br), bindex_(be) {}

        const Molecule& br() const {
            return *br_;
        }

        bool operator==(const BondIterator& rhs) const;

        bool operator!=(const BondIterator& rhs) const;

        BondIterator& operator++();

        BondIterator operator++(int);

        BondEdge operator* () const;
    private:
        const Molecule* br_;
        OutEdgeIterator bindex_;
    };

    Bonds(OutBondIteratorPair be, const Molecule* br) :
        begin_end_(be), br_(br) {}

    BondIterator begin() const {
        return {begin_end_.first, br_};
    }

    BondIterator end() const {
        return {begin_end_.second, br_};
    }
private:
    const OutBondIteratorPair begin_end_;
    const Molecule* br_;
};

class SPEAR_EXPORT BondEdge {
public:
    BondEdge(const Molecule* br, EdgeDescriptor index) : br_(br), index_(index) {}

    AtomVertex source() const;

    AtomVertex target() const;

    AtomVertex other_atom(const AtomVertex& av) const;

    std::pair<AtomVertex, AtomVertex> as_pair() const;

    size_t index() const;

    Bond::Order order() const;

private:
    const Molecule* br_;
    EdgeDescriptor index_;
};

class SPEAR_EXPORT AtomVertex {
public:

    AtomVertex(const Molecule* br, VertexDescriptor index);

    const Molecule& br() const {
        return *br_;
    }

    const std::string& name() const;
    
    const Eigen::Vector3d& position() const;

    Element::Symbol atomic_number() const;

    size_t degree() const;
    
    Neighbors neighbors() const;

    Bonds bonds() const;

    bool is_aromatic() const;

    bool is_planar() const;

    bool is_non_metal() const;

    size_t expected_bonds() const;

    size_t implicit_hydrogens() const;

    size_t explicit_hydrogens() const;

    size_t total_hydrogens() const;

    Atom::Chirality chirality() const;

    AtomRingMapIteratorPair rings() const;

    AtomRingMapIteratorPair sssrs() const;

    Hybridization hybridization() const;

    double partial_charge() const;

    int8_t formal_charge() const; 

    bool operator==(const AtomVertex& rhs) const;

    AtomVertex operator[](size_t) const;

    operator size_t() const;
    
private:
    VertexDescriptor index_;
    const Molecule* br_;
};

}

#include "spear/Graph_impl.hpp"

#endif
