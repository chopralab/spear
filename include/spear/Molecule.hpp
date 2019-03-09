// License:

#ifndef SPEAR_MOLECULE_HPP
#define SPEAR_MOLECULE_HPP

#include "spear/exports.hpp"

#include <vector>

#include "chemfiles/Frame.hpp"
#include "chemfiles/external/optional.hpp"

#include "spear/Graph.hpp"
#include "spear/AtomType.hpp"

namespace Spear {

class SPEAR_EXPORT Molecule {
public:

    class iterator {
    public:

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::random_access_iterator_tag;
        using value_type = AtomVertex;
        using pointer = AtomVertex*;
        using reference = AtomVertex&;

        iterator() : index_(0), base_mol_(nullptr) {}

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

    explicit Molecule(chemfiles::Frame frame) :
        frame_(std::move(frame)), graph_() {
        init_();
    }

    const Graph& graph() const {
        return graph_;
    }

    const chemfiles::Frame& frame() const {
        return frame_;
    }

    const std::set<std::set<size_t>> rings() const;

    size_t dimensionality(double eps = 0.00001) const;

    std::vector<EdgeDescriptor> get_bonds_in(const std::set<size_t>& atoms) const;

    template<class atomtype, typename... args>
    std::string add_atomtype(args... additional);

    optional<const AtomType*> get_atomtype(const std::string& name) const;

    void set_default_atomtype(const std::string& name);

    const AtomType* get_default_atomtype() const;

    void remove_hydrogens();

    size_t size() const;

    AtomVertex operator[](size_t index) const;

    iterator begin() const;

    iterator end() const;

    iterator cbegin() const;

    iterator cend() const;

private:
    void init_();

    chemfiles::Frame frame_;
    
    Graph graph_;

    typedef std::unique_ptr<AtomType> unqiue_AtomType;
    std::unordered_map<std::string, unqiue_AtomType> atom_types_;

    std::string default_atomtype_;
};

}

#include "spear/Molecule_impl.hpp"

#endif
