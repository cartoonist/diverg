/**
 *    @file  adjacency.hpp
 *   @brief  Graph adjacency matrix construction (Kokkos-free).
 *
 *  This header offers interface functions to build the (Basic/Range) adjacency
 *  matrix of a sequence graph, plus the small graph helpers it needs.
 *
 *  This is *host-only* code that produces a Kokkos-free `CRSMatrix`, so it can
 *  be included by downstream host translation units without pulling in Kokkos.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_ADJACENCY_HPP__
#define DIVERG_ADJACENCY_HPP__

#include <cassert>
#include <cstddef>
#include <algorithm>
#include <string>
#include <tuple>

#include <gum/graph.hpp>

#include "crs_matrix.hpp"


namespace diverg {
  namespace util {
    template< typename TGraph >
    inline auto
    first_node_in_region( TGraph const& graph,
                          std::string const& reg_name )
    {
      typename TGraph::id_type first = 0;

      graph.for_each_path( [&]( auto, auto id ) {
        if ( graph.path_name( id ) == reg_name ) {
          auto path = graph.path( id );
          first = path.id_of( path.front() );
          return false;
        }
        return true;
      } );

      return first;
    }

    template< typename TGraph >
    inline auto
    node_by_name( TGraph const& graph, std::string const& name )
    {
      typename TGraph::id_type nodeid = 0;

      graph.for_each_node( [&]( auto rank, auto id ) {
        if ( graph.get_node_prop( rank ).name == name ) {
          nodeid = id;
          return false;
        }
        return true;
      } );

      return nodeid;
    }

    /**
     *  @brief  Get the upper bound on node ranks in a component containing the
     *          given node with rank `node_rank`.
     *
     *  @param[in]  graph      The graph.
     *  @param[in]  node_rank  The rank of the node.
     *
     *  @return  A tuple containing the upper bound on node ranks in the
     *           component, component index the `node_rank` belongs to, and the
     *           total number of components.
     *
     *  NOTE: This method assumes that the input graph is sorted such that node
     *  rank ranges in components are disjoint.
     */
    template< typename TGraph >
    inline auto
    upper_node_rank_in_component( TGraph const& graph,
                                  typename TGraph::rank_type node_rank )
    {
      using rank_type = typename TGraph::rank_type;

      rank_type upper_rank = 0;
      std::size_t nof_comps = 0;
      rank_type comp_i = 0;
      gum::util::for_each_start_side( graph, [&]( auto rank, auto ) {
        ++nof_comps;
        if ( upper_rank == 0 ) {
          if ( rank > node_rank )
            upper_rank = rank;
          else
            comp_i = nof_comps;
        }
        return true;
      } );

      if ( upper_rank == 0 ) upper_rank = graph.get_node_count() + 1;

      return std::make_tuple( upper_rank, comp_i, nof_comps );
    }

    /**
     *  @brief  Count the nodes of a subgraph in the graph or of the whole graph.
     *
     *  The subgraph is indicated by the node rank range [lower, upper).
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers
     *  all nodes in the component.
     */
    template< class TGraph >
    inline typename TGraph::rank_type
    _node_count( TGraph const& graph, typename TGraph::rank_type lower=1,
                typename TGraph::rank_type upper=0 )
    {
      if ( upper == 0 ) upper = graph.get_node_count() + 1;
      assert( lower <= graph.get_node_count() && lower > 0 );
      assert( upper <= graph.get_node_count()+1 && upper > lower );
      return upper - lower;
    }

    /**
     *  @brief  Count the edges of a graph component or of the whole graph.
     *
     *  The component is indicated by the node range [lower, upper).
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers
     *  all nodes in the component.
     */
    template< class TGraph >
    inline typename TGraph::rank_type
    _edge_count( TGraph const& graph, typename TGraph::rank_type lower=1,
                typename TGraph::rank_type upper=0 )
    {
      typedef typename TGraph::id_type id_type;
      typedef typename TGraph::rank_type rank_type;

      rank_type count = 0;
      graph.for_each_node(
          [&graph, &count, &upper]( rank_type rank, id_type id ){
            count += graph.outdegree( id );
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      return count;
    }

    /**
     *  @brief  Get the adjacency matrix of the graph as a native Basic CRS matrix.
     *
     *  @param[in]  graph The graph.
     *  @param[in]  lower  The lower node rank (inclusive).
     *  @param[in]  upper  The upper node rank (exclusive).
     *  @return  The adjacency matrix as a native Basic `TCRSMatrix`.
     *
     *  Overload of `adjacency_matrix` for native `CRSMatrix` types. It builds
     *  the matrix for the subgraph defined by the node rank range [lower,
     *  upper). It uses Kokkos core only, decoupling from KokkosKernels.
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers
     *  all nodes in the component.
     */
    template< typename TCRSMatrix, class TGraph,
              std::enable_if_t<
                  crs_matrix::is_native_crs_matrix< TCRSMatrix >::value, int > = 0 >
    inline TCRSMatrix
    adjacency_matrix( TGraph const& graph,
                      typename TGraph::rank_type lower=1,
                      typename TGraph::rank_type upper=0 )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::offset_type offset_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::linktype_type linktype_type;

      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;
      typedef typename crsmat_type::size_type size_type;

      static_assert( crs_matrix::is_basic_crs_matrix< crsmat_type >::value,
                     "output matrix should be in Basic CRS format." );

      if ( upper == 0 ) upper = graph.get_node_count() + 1;
      ordinal_type nrows = gum::util::total_nof_loci( graph, lower, upper );
      size_type nnz = nrows - _node_count( graph, lower, upper ) +
          _edge_count( graph, lower, upper );

      auto entries = crsmat_type::make_entries( nnz );
      auto rowmap = crsmat_type::make_rowmap( nrows + 1 );

      offset_type cursor = 0;
      offset_type start = gum::util::id_to_charorder( graph, graph.rank_to_id( lower ) );
      size_type i = 0;
      size_type irow = 0;
      rowmap[ irow++ ] = i;
      graph.for_each_node(
          [&]( rank_type rank, id_type id ) {
            for ( offset_type offset = 1; offset < graph.node_length( id ); ++offset ) {
              entries[ i++ ] = ++cursor;
              rowmap[ irow++ ] = i;
            }
            ++cursor;
            auto row_begin = i;
            graph.for_each_edges_out(
                id,
                [&graph, &entries, &i, start]( id_type to, linktype_type ) {
                  entries[ i++ ] = gum::util::id_to_charorder( graph, to ) - start;
                  return true;
                } );
            std::sort( entries.data() + row_begin, entries.data() + i );
            rowmap[ irow++ ] = i;
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      assert( i == nnz );
      assert( irow == static_cast< size_type >( nrows + 1 ) );

      return crsmat_type( nrows, std::move( entries ), std::move( rowmap ) );
    }

    /**
     *  @brief  Get the adjacency matrix of the graph in Range CRS format.
     *
     *  @param[in]  graph  The graph.
     *  @param[in]  lower  The lower node rank (inclusive).
     *  @param[in]  upper  The upper node rank (exclusive).
     *  @return  The adjacency matrix as a Range CRS matrix (`TRangeCRSMatrix`).
     *
     *  Like `adjacency_matrix`, but yields the result in Range CRS format using
     *  only Kokkos core; it does not go through a `KokkosSparse::CrsMatrix` and
     *  therefore does not require KokkosKernels.
     *
     *  It builds the adjacency in the corresponding native Basic format (via the
     *  native `adjacency_matrix` overload) and converts it to Range format by
     *  assignment, which performs the basic-to-range conversion internally.
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers all
     *  nodes in the component(s).
     */
    template< typename TRangeCRSMatrix, class TGraph >
    inline TRangeCRSMatrix
    range_adjacency_matrix( TGraph const& graph,
                            typename TGraph::rank_type lower=1,
                            typename TGraph::rank_type upper=0 )
    {
      typedef TRangeCRSMatrix crsmat_type;

      static_assert( crs_matrix::is_range_crs_matrix< crsmat_type >::value,
                     "output matrix should be in Range CRS format." );

      typedef diverg::make_basic_t< crsmat_type > basic_crsmat_type;

      crsmat_type matrix;
      matrix.assign( adjacency_matrix< basic_crsmat_type >( graph, lower, upper ) );
      return matrix;
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_ADJACENCY_HPP__ --- */
