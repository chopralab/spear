// Spear: Statistical Platform for Elucidating moleculAr Reactivity
// Copyright (C) Purdue University -- BSD license

#ifndef SPEAR_TYPEDEFS_HPP
#define SPEAR_TYPEDEFS_HPP

#include <boost/graph/undirected_graph.hpp>
#include "spear/Constants.hpp"

namespace Spear {

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

}

#endif
