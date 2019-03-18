#ifndef SPEAR_GRAPH_HPP
#define SPEAR_GRAPH_HPP

#include <boost/graph/undirected_graph.hpp>

#include "chemfiles/external/optional.hpp"

#include "spear/Constants.hpp"
#include "spear/Rings.hpp"

namespace Spear {

using chemfiles::optional;

// Note: If more atom/bond properties are need other than atomic number and bond order,
// change the lines below:
using VertexProperty = boost::property<boost::vertex_name_t, Element::Symbol>;
using EdgeProperty = boost::property<boost::edge_name_t, Bond::Order>;

// setS is used because we do **NOT** want parallel graphs (bonds represented multiple times)
// vecS is used because we do **WANT** verticies to be one-to-one with their index
// Molecules are naturally undirected, with the exception of dative bonds (handled by the edge property)
using Graph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS,
                                    VertexProperty, EdgeProperty>;

using Traits = boost::graph_traits<Graph>;
using VertexDescriptor = Traits::vertex_descriptor;
using EdgeDescriptor = Traits::edge_descriptor;

using VertexIterator = Traits::vertex_iterator;
using AdjacencyIterator = Traits::adjacency_iterator;
using EdgeIterator = Traits::edge_iterator;
using OutEdgeIterator = Traits::out_edge_iterator;

using AdjacencyIteratorPair = std::pair<AdjacencyIterator, AdjacencyIterator>;
using BondIteratorPair = std::pair<EdgeIterator, EdgeIterator>;
using OutBondIteratorPair = std::pair<OutEdgeIterator, OutEdgeIterator>;

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
            bindex_(be), br_(br) {}

        bool operator==(const BondIterator& rhs) const;

        bool operator!=(const BondIterator& rhs) const;

        BondIterator& operator++();

        BondIterator operator++(int);

        BondEdge operator* () const;
    private:
        OutEdgeIterator bindex_;
        const Molecule* br_;
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

    Bond::Order order() const;

private:
    const Molecule* br_;
    EdgeDescriptor index_;
};

class SPEAR_EXPORT AtomVertex {
public:

    AtomVertex(const Molecule* br, VertexDescriptor index);

    const std::string& name() const;

    const std::string& type() const;
    
    const Eigen::Vector3d& position() const;

    Element::Symbol atomic_number() const;

    size_t degree() const;
    
    Neighbors neighbors() const;

    Bonds bonds() const;

    bool is_aromatic() const;

    bool is_non_metal() const;

    size_t expected_bonds() const;

    size_t implicit_hydrogens() const;

    size_t explicit_hydrogens() const;

    size_t total_hydrogens() const;

    AtomRingMapIteratorPair rings() const;

    AtomRingMapIteratorPair sssrs() const;

    bool operator==(const AtomVertex& rhs) const;

    AtomVertex operator[](size_t) const;

    operator size_t() const;
    
private:
    VertexDescriptor index_;
    const Molecule* br_;
};

}

#endif
