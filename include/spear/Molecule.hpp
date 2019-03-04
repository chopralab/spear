// License:

#ifndef SPEAR_MOLECULE_HPP
#define SPEAR_MOLECULE_HPP

#include "spear/exports.hpp"

#include <vector>

#include "chemfiles/Frame.hpp"
#include "chemfiles/external/optional.hpp"

#include <boost/graph/undirected_graph.hpp>

#include "spear/AtomType.hpp"

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

using chemfiles::optional;

class AtomVertex {
public:
    AtomVertex(const Molecule* br, size_t index);

    const std::string& name() const;

    const std::string& type() const;
    
    const chemfiles::Vector3D& position() const;

    size_t atomic_number() const;

    size_t neighbor_count() const;
    
    AdjacencyIteratorPair neighbors() const;

    bool is_aromatic() const;

    size_t expected_bonds() const;

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

    class iterator {
    public:

        using difference_type = std::ptrdiff_t;

        using iterator_category = std::random_access_iterator_tag;
        using value_type = AtomVertex;
        using pointer = AtomVertex*;
        using reference = AtomVertex&;

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

    std::vector<EdgeDescriptor> get_bonds_in(const std::set<size_t>& atoms) const;

    template<class atomtype, typename typemode>
    std::string add_atomtype(typemode mode);

    optional<const AtomType*> get_atomtype(const std::string& name) const;

    void set_default_atomtype(const std::string& name);

    const AtomType* get_default_atomtype() const;

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
    return static_cast<size_t>(std::distance(tmp.first, tmp.second));
}

inline AdjacencyIteratorPair AtomVertex::neighbors() const {
    return boost::adjacent_vertices(index_, br_->graph());
}

inline bool AtomVertex::is_aromatic() const {
    auto types = br_->get_default_atomtype();
    return types->is_aromatic(index_);
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

template<class atomtype, typename typemode>
inline std::string Molecule::add_atomtype(typemode mode) {
    auto typed_atoms = new atomtype(*this, mode);
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

    if (default_atomtype_.empty()) {
        throw std::runtime_error("No atom types are availible.");
    }

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
