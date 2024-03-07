/**
 *    @file  dindex.hpp
 *   @brief  Distance index module.
 *
 *  This header file contains all APIs, general utility and helper functions
 *  related to distance index.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Wed Mar 06, 2024  18:56
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2024, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_DINDEX_HPP__
#define DIVERG_DINDEX_HPP__

#include <algorithm>

#include <gum/graph.hpp>


namespace diverg {
  namespace util {
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
     *  @brief  Get the adjacency matrix of the graph in CRS format.
     *
     *  @param[in]  graph The graph.
     *  @param[in]  lower  The lower node rank (inclusive).
     *  @param[in]  upper  The upper node rank (exclusive).
     *  @return  The adjacency matrix in KokkosSparse::CrsMatrix format.
     *
     *  Compute adjacency matrix of a component in the given `graph` or of the whole
     *  graph. The component is indicated by nodes whose ranks are in the range [lower,
     *  upper). The resulting adjacency matrix is stored in CRS format.
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers
     *  all nodes in the component.
     */
    template< typename TCrsMatrix, class TGraph >
    inline TCrsMatrix
    adjacency_matrix( TGraph const& graph,
                      typename TGraph::rank_type lower=1,
                      typename TGraph::rank_type upper=0 )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::id_type id_type;
      typedef typename graph_type::offset_type offset_type;
      typedef typename graph_type::rank_type rank_type;
      typedef typename graph_type::linktype_type linktype_type;

      typedef typename TCrsMatrix::staticcrsgraph_type staticcrsgraph_type;
      typedef typename TCrsMatrix::values_type::non_const_type values_type;
      typedef typename staticcrsgraph_type::data_type ordinal_type;
      typedef typename staticcrsgraph_type::size_type size_type;
      typedef typename staticcrsgraph_type::row_map_type::non_const_type row_map_type;
      typedef typename staticcrsgraph_type::entries_type::non_const_type entries_type;

      if ( upper == 0 ) upper = graph.get_node_count() + 1;
      ordinal_type nrows = gum::util::total_nof_loci( graph, lower, upper );
      size_type nnz = nrows - _node_count( graph, lower, upper ) +
          _edge_count( graph, lower, upper );

      entries_type entries( "entries", nnz );
      values_type values( "values", nnz );
      row_map_type rowmap( "rowmap", nrows + 1 );

      for ( size_type i = 0; i < nnz; ++i ) values( i ) = 1;  // boolean

      offset_type cursor = 0;
      offset_type start = gum::util::id_to_charorder( graph, graph.rank_to_id( lower ) );
      size_type i = 0;
      size_type irow = 0;
      rowmap( irow++ ) = i;
      graph.for_each_node(
          [&]( rank_type rank, id_type id ) {
            assert( gum::util::id_to_charorder( graph, id ) == cursor + start );
            for ( offset_type offset = 1; offset < graph.node_length( id ); ++offset ) {
              entries( i++ ) = ++cursor;
              rowmap( irow++ ) = i;
            }
            ++cursor;
            auto entries_begin = entries.data() + i;
            graph.for_each_edges_out(
                id,
                [&graph, &entries, &i, start]( id_type to, linktype_type ) {
                  entries( i++ ) = gum::util::id_to_charorder( graph, to ) - start;
                  return true;
                } );
            std::sort( entries_begin, entries.data() + i );
            rowmap( irow++ ) = i;
            if ( rank + 1 == upper ) return false;
            return true;
          },
          lower );
      assert( i == nnz );
      assert( irow == static_cast< unsigned int >( nrows + 1 ) );

      return TCrsMatrix( "adjacency matrix", nrows, nrows, nnz, values, rowmap, entries );
    }

    /**
     *  @brief  Compress a distance index by removing intra-node loci pairs
     *
     *  @param  dindex input distance index
     *  @param  graph underlying graph
     *  @return a mutable compressed distance index of type `TMutableCRSMatrix`
     *
     *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
     *        matrix afterwards.
     *
     *  NOTE: The input uncompressed distance index is passed by non-const reference,
     *        since containers in const Buffered specialisations cannot be iterated.
     */
    template< typename TMutableCRSMatrix,
              typename TCRSMatrix,
              typename TGraph >
    inline TMutableCRSMatrix
    compress_distance_index( TCRSMatrix& dindex, TGraph const& graph )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::rank_type rank_type;

      typedef TMutableCRSMatrix crsmat_mutable_type;
      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;
      typedef typename crsmat_type::size_type size_type;

      typename crsmat_mutable_type::entries_type entries;
      typename crsmat_mutable_type::rowmap_type rowmap;
      crsmat_mutable_type::base_type::traits_type::init( entries );
      crsmat_mutable_type::base_type::traits_type::init( rowmap );
      rank_type cnode_rank = 0;  // current node rank
      size_type start = 0;    // row start index
      size_type end;          // row end index
      ordinal_type cloc = 0;     // current node loci index
      ordinal_type nloc = 0;     // next node loci index
      for ( ordinal_type nrow = 0; nrow < dindex.numRows(); ++nrow ) {
        rowmap.push_back( entries.size() );
        if ( nrow == nloc ) {
          ++cnode_rank;
          cloc = nloc;
          nloc += graph.node_length( graph.rank_to_id( cnode_rank ) );
        }
        assert( nrow < nloc );
        end = dindex.rowMap( nrow + 1 );
        for ( ; start < end; ++start ) {
          if ( cloc <= dindex.entry( start ) && dindex.entry( start ) < nloc ) continue;
          else entries.push_back( dindex.entry( start ) );
        }
      }
      rowmap.push_back( entries.size() );
      assert( start == dindex.nnz() );

      return crsmat_mutable_type( dindex.numCols(), std::move( entries ), std::move( rowmap ) );
    }
  } // namespace util
} // namespace diverg

#endif // DIVERG_DINDEX_HPP__
