// License:

#ifndef SPEAR_MOLECULE_HPP
#define SPEAR_MOLECULE_HPP

#include "spear/exports.hpp"

#include <vector>

#include "chemfiles/Frame.hpp"

#include "Eigen/Geometry"

#include "spear/Graph.hpp"
#include "spear/Rings.hpp"
#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT Molecule {
public:

    template<typename ret_type, typename sto_type>
    class iterator {
    public:

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::random_access_iterator_tag;
        using value_type = ret_type;
        using pointer = ret_type*;
        using reference = ret_type&;

        iterator() : index_(0), base_mol_(nullptr) {}

        iterator(const Molecule* base_mol, sto_type i) : index_(i), base_mol_(base_mol) {
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

        value_type operator* () const;

        value_type operator[](difference_type rhs) const;

    private:
        sto_type index_;
        const Molecule* base_mol_;
    };

    struct AllBonds {
        AllBonds(const Molecule* br, BondIteratorPair be) :
            begin_end_(be), br_(br) {}

        iterator<BondEdge, EdgeIterator> begin() const {
            return iterator<BondEdge, EdgeIterator>(br_, begin_end_.first);
        }

        iterator<BondEdge, EdgeIterator> end() const {
            return iterator<BondEdge, EdgeIterator>(br_, begin_end_.second);
        }

        iterator<BondEdge, EdgeIterator> cbegin() const {
            return iterator<BondEdge, EdgeIterator>(br_, begin_end_.first);
        }

        iterator<BondEdge, EdgeIterator> cend() const {
            return iterator<BondEdge, EdgeIterator>(br_, begin_end_.second);
        }

    private:
        BondIteratorPair begin_end_;
        const Molecule* br_;
    };

    explicit Molecule(chemfiles::Frame frame) :
        frame_(std::move(frame)), graph_() {
        init_();
        rings_();
        smallest_set_of_smallest_rings_();
    }

    const Graph& graph() const {
        return graph_;
    }

    const chemfiles::Frame& frame() const {
        return frame_;
    }

    const std::vector<Eigen::Vector3d>& positions() const {
        return positions_;
    }

    const RingSet& rings() const {
        return all_rings_;
    }

    const RingSet& smallest_set_of_smallest_rings() const {
        return sssr_;
    }

    size_t dimensionality(double eps = 0.00001) const;

    BondEdge bond(size_t idx1, size_t idx2);

    std::vector<BondEdge> get_bonds_in(const std::set<size_t>& atoms) const;

    template<class atomtype, typename... args>
    std::string add_atomtype(args... additional);

    const AtomType* get_atomtype(const std::string& name) const;

    void set_default_atomtype(const std::string& name);

    const AtomType* get_default_atomtype() const;

    AtomVertex add_atom(Element::Symbol n_atom, const Eigen::Vector3d& pos);

    BondEdge add_bond(size_t idx1, size_t idx2, Bond::Order order = Bond::SINGLE);

    AtomVertex add_atom_to(Element::Symbol n_atom, size_t index);

    void remove_hydrogens();

    size_t add_hydrogens();

    size_t size() const;

    AtomVertex operator[](size_t index) const;

    iterator<AtomVertex, VertexIterator> begin() const;

    iterator<AtomVertex, VertexIterator> end() const;

    iterator<AtomVertex, VertexIterator> cbegin() const;

    iterator<AtomVertex, VertexIterator> cend() const;

    AllBonds bonds() const;

private:
    void init_();

    void rings_();

    void smallest_set_of_smallest_rings_();

    /// Chemfiles frame
    chemfiles::Frame frame_;

    /// Coordinates of the molecule
    std::vector<Eigen::Vector3d> positions_;

    /// Graph representation of the molecule
    Graph graph_;

    /// All rings found in the molecule
    RingSet all_rings_;

    /// Smallest set of smallest rings
    RingSet sssr_;

    /// Maps an atom to all rings that contain it
    AtomRingMap atom_to_ring_;

    /// Maps an atom to all sssr rings that contain it
    AtomRingMap atom_to_sssr_;

    typedef std::unique_ptr<AtomType> unqiue_AtomType;
    std::unordered_map<std::string, unqiue_AtomType> atom_types_;

    std::string default_atomtype_;

    friend class AtomVertex;
};

}

#endif
