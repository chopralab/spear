// License:

#ifndef SPEAR_MOLECULE_HPP
#define SPEAR_MOLECULE_HPP

#include "spear/exports.hpp"

#include <vector>

#include "chemfiles/Frame.hpp"

#include "Eigen/Geometry"

#include "spear/Typedefs.hpp"
#include "spear/Rings.hpp"
#include "spear/AtomType.hpp"
#include "spear/PartialCharge.hpp"

namespace Spear {

class BondEdge;
class AllBonds;
class AtomVertex;

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

    Molecule();

    explicit Molecule(const chemfiles::Frame& frame) :
        graph_(), topology_(frame.topology()) {
        init_(frame);
        rings_();
        add_default_type();
    }

    const chemfiles::Topology& topology() const {
        return topology_;
    }

    const Graph& graph() const {
        return graph_;
    }

    const std::vector<Eigen::Vector3d>& positions() const {
        return positions_;
    }

    void set_positions(std::vector<Eigen::Vector3d> pos) {
        positions_ = std::move(pos);
    }

    const RingSet& rings() const {
        return all_rings_;
    }

    const RingSet& smallest_set_of_smallest_rings() const {
        return sssr_;
    }

    /// Retreive a bond in the molecule given the two indicies.
    ///
    /// \param idx1 first atom index
    /// \param idx2 second atom index
    BondEdge bond(size_t idx1, size_t idx2) const;

    /// Retrieve all the bonds within all the atoms.
    ///
    /// \param atoms A set of atom indicies
    /// \returns A vector of `BondEdge` which are always between the `atoms`.
    std::vector<BondEdge> get_bonds_in(const std::set<size_t>& atoms) const;

    /// Determine the number of connected graphs in the molecule
    ///
    /// \returns the number of independant graphs
    size_t connected_components() const;

    /// Add an `AtomType` to the `Molecule`.
    ///
    /// \params additional The arguments to constructor of `atomtype`.
    /// \returns The name of the new `AtomType`.
    template<class atomtype, typename... args>
    std::string add_atomtype(args... additional);

    /// Retreive an `AtomType` for the `Molecule`
    ///
    /// \params name The name of the `AtomType` to retrieve.
    /// \returns A pointer to the `AtomType` or nullptr if the name is not found
    const AtomType* atomtype(const std::string& name = "") const;

    /// Set the default `AtomType` 
    ///
    /// The definitions of aromaticity, hyridization, and charge are dependant
    /// on the current `AtomType` beging used. While molecules can have
    /// multiple `AtomType` assignments, they can only have one default
    /// assignment. This default is used for graph based searches where these
    /// properties are used. This function changes this default atom type.
    /// \params name The name of the `AtomType` to be used as the default.
    void set_default_atomtype(const std::string& name);

    /// Assign a `PartialCharge` to the `Molecule`.
    ///
    /// \params additional The arguments to constructor of `partialcharge`.
    /// \returns The name of the new `PartialCharge`.
    template<class partialcharge, typename... args>
    std::string add_partial_charge(args... additional);

    /// Retreive a `PartialCharge` for the `Molecule`
    ///
    /// \params name The name of the `PartialCharge` to retrieve.
    /// \returns A pointer to the `PartialCharge` or nullptr if the name is not found
    const PartialCharge* partial_charge(const std::string& name = "") const;

    /// Set the default `PartialCharge`
    ///
    /// The definitions of aromaticity, hyridization, and charge are dependant
    /// on the current `PartialCharge` beging used. While molecules can have
    /// multiple `PartialCharge` assignments, they can only have one default
    /// assignment. This default is used for graph based searches where these
    /// properties are used. This function changes this default atom type.
    /// \params name The name of the `PartialCharge` to be used as the default.
    void set_default_partial_charge(const std::string& name);

    AtomVertex add_atom(Element::Symbol n_atom, const Eigen::Vector3d& pos);

    BondEdge add_bond(size_t idx1, size_t idx2, Bond::Order order = Bond::SINGLE);

    AtomVertex add_atom_to(Element::Symbol n_atom, size_t index);

    void remove_bond(size_t idx1, size_t idx2);

    void remove_atom(size_t index);

    void swap_atoms(size_t idx1, size_t idx2);

    void remove_hydrogens();

    size_t add_hydrogens();

    /// The number of atoms(nodes) in the `Molecule`
    size_t size() const;

    /// The number of bonds(edges) in the `Molecule`
    size_t bond_count() const;

    AtomVertex operator[](size_t index) const;

    iterator<AtomVertex, VertexIterator> begin() const;

    iterator<AtomVertex, VertexIterator> end() const;

    iterator<AtomVertex, VertexIterator> cbegin() const;

    iterator<AtomVertex, VertexIterator> cend() const;

    AllBonds bonds() const;

private:
    /// Initialize a `Molecule` from a chemfiles frame
    void init_(const chemfiles::Frame& frame);

    /// Calculate the rings of a molecule
    void rings_();

    /// Calculate the smallest set of smallest rings
    void smallest_set_of_smallest_rings_();

    /// Adds the default type after construction
    void add_default_type();

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

    /// Store just the topology
    chemfiles::Topology topology_;

    /// Storage of all `AtomType` assigned for the molecule
    typedef std::unique_ptr<AtomType> unqiue_AtomType;
    std::unordered_map<std::string, unqiue_AtomType> atom_types_;

    /// The current default `AtomType` 
    std::string default_atomtype_;

    /// Storage of all PartialCharge assigned for the molecule
    typedef std::unique_ptr<PartialCharge> unique_PartialCharge;
    std::unordered_map<std::string, unique_PartialCharge> partial_charges_;

    /// The current default `PartialCharge`
    std::string default_partial_charge_;

    friend class AtomVertex;
};

}

#include "spear/Molecule_impl.hpp"

#endif
