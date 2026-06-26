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
#include <vector>

#include <gum/graph.hpp>

#include "range_sparse.hpp"


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
     *  @brief  Get the adjacency matrix of the graph as a non-native CRS matrix.
     *
     *  @param[in]  graph The graph.
     *  @param[in]  lower  The lower node rank (inclusive).
     *  @param[in]  upper  The upper node rank (exclusive).
     *  @return  The adjacency matrix as a `KokkosSparse::CrsMatrix`-like matrix.
     *
     *  Compute adjacency matrix of a component in the given `graph` or of the whole
     *  graph. The component is indicated by nodes whose ranks are in the range [lower,
     *  upper). The resulting adjacency matrix is stored in CRS format.
     *
     *  This overload targets non-native output types (e.g. `KokkosSparse::CrsMatrix`);
     *  see the native overload for `CRSMatrix` types.
     *
     *  NOTE: This function assumes the range [lower, upper) exclusively covers
     *  all nodes in the component.
     */
    template< typename TCrsMatrix, class TGraph,
              std::enable_if_t< !crs_matrix::is_native_crs_matrix< TCrsMatrix >::value, int > = 0 >
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
      assert( irow == static_cast< size_type >( nrows + 1 ) );

      return TCrsMatrix( "adjacency matrix", nrows, nrows, nnz, values, rowmap, entries );
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

    /**
     *  @brief  Merge two distance indices (Range Group, out-of-core operands)
     *
     *  @param  dindex1 first distance index
     *  @param  dindex2 second distance index
     *  @return a mutable merged distance index of type `TMutableCRSMatrix`
     *
     *  Host-sequential streaming merge used when an input or the output is
     *  buffered (out-of-core), where `range_spadd` (which requires full
     *  residency) cannot be applied.
     *
     *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
     *        matrix afterwards.
     *
     *  NOTE: The input distance indices are passed by non-const references, since
     *        containers in const Buffered specialisations cannot be iterated.
     */
    template< typename TMutableCRSMatrix, typename TCRSMatrix >
    inline TMutableCRSMatrix
    merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                          crs_matrix::RangeGroup /* group tag */,
                          std::true_type /* is_buffered */ )
    {
      typedef TMutableCRSMatrix crsmat_mutable_type;
      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;
      typedef typename crsmat_type::size_type size_type;

      auto entries = crsmat_mutable_type::make_entries();
      auto rowmap = crsmat_mutable_type::make_rowmap();
      size_type cursor1 = 0;    // current entry index in the first distance index
      size_type cursor2 = 0;    // current entry index in the second distance index
      size_type end1;  // last entry index of the row in the first distance index
      size_type end2;  // last entry index of the row in the second distance index

      assert( dindex1.numRows() == dindex2.numRows() );
      assert( dindex1.numCols() == dindex2.numCols() );

      auto nof_rows = dindex1.numRows();
      auto nof_cols = dindex1.numCols();

      /* NOTE: this function assumes that `cursors` are in the entry ranges. */
      auto fetch_min_and_adv =
        []( TCRSMatrix& dindex1, TCRSMatrix& dindex2, size_type& cursor1, size_type& cursor2 )
        {
          ordinal_type l1 = dindex1.entry( cursor1 );
          ordinal_type u1 = dindex1.entry( cursor1 + 1 );
          ordinal_type l2 = dindex2.entry( cursor2 );
          ordinal_type u2 = dindex2.entry( cursor2 + 1 );

          if ( l1 < l2 || ( l1 == l2 && u1 < u2 ) ) {
            ++cursor1; ++cursor1;
            return std::make_pair( l1, u1 );
          }
          else {
            ++cursor2; ++cursor2;
            return std::make_pair( l2, u2 );
          }
        };

      auto merge_and_adv =
        []( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
            size_type& cursor1, size_type end1,
            size_type& cursor2, size_type end2,
            ordinal_type l, ordinal_type u )
        {
          ordinal_type l1 = std::numeric_limits< ordinal_type >::max();
          ordinal_type u1 = std::numeric_limits< ordinal_type >::max();
          ordinal_type l2 = std::numeric_limits< ordinal_type >::max();
          ordinal_type u2 = std::numeric_limits< ordinal_type >::max();

          if ( cursor1 < end1 ) {
            l1 = dindex1.entry( cursor1 );
            u1 = dindex1.entry( cursor1 + 1 );
          }

          if ( cursor2 < end2 ) {
            l2 = dindex2.entry( cursor2 );
            u2 = dindex2.entry( cursor2 + 1 );
          }

          while ( true ) {
            if ( u + 1 >= l1 ) {
              l = std::min( l, l1 );
              u = std::max( u, u1 );
              ++cursor1; ++cursor1;
              if ( cursor1 < end1 ) {
                l1 = dindex1.entry( cursor1 );
                u1 = dindex1.entry( cursor1 + 1 );
              }
              else l1 = std::numeric_limits< ordinal_type >::max();
            }
            else if ( u + 1 >= l2 ) {
              l = std::min( l, l2 );
              u = std::max( u, u2 );
              ++cursor2; ++cursor2;
              if ( cursor2 < end2 ) {
                l2 = dindex2.entry( cursor2 );
                u2 = dindex2.entry( cursor2 + 1 );
              }
              else l2 = std::numeric_limits< ordinal_type >::max();
            }
            else break;
          }

          return std::make_pair( l, u );
        };

      size_type c_nnz = 0;
      ordinal_type l = 0;
      ordinal_type u = 0;
      for ( ordinal_type nrow = 0; nrow < nof_rows; ++nrow ) {
        rowmap.push_back( entries.size() );

        end1 = dindex1.rowMap( nrow + 1 );
        end2 = dindex2.rowMap( nrow + 1 );
        l = 0;
        u = 0;
        while ( cursor1 < end1 ) {
          if ( cursor2 >= end2 ) {
            while ( cursor1 < end1 ) {
              l = dindex1.entry( cursor1++ );
              u = dindex1.entry( cursor1++ );
              entries.push_back( l );
              entries.push_back( u );
              c_nnz += u - l + 1;
            }
            break;
          }
          std::tie( l, u ) = fetch_min_and_adv( dindex1, dindex2, cursor1, cursor2 );
          std::tie( l, u ) =
            merge_and_adv( dindex1, dindex2, cursor1, end1, cursor2, end2, l, u );
          entries.push_back( l );
          entries.push_back( u );
          c_nnz += u - l + 1;
        }
        while ( cursor2 < end2 ) {
          l = dindex2.entry( cursor2++ );
          u = dindex2.entry( cursor2++ );
          entries.push_back( l );
          entries.push_back( u );
          c_nnz += u - l + 1;
        }
      }
      rowmap.push_back( entries.size() );

      return crsmat_mutable_type( nof_cols, std::move( entries ), std::move( rowmap ), c_nnz );
    }

    /**
     *  @brief  Merge two distance indices (Range Group, in-core operands)
     *
     *  @param  dindex1 first distance index
     *  @param  dindex2 second distance index
     *  @return a mutable merged distance index of type `TMutableCRSMatrix`
     *
     *  Merging two Boolean distance indices is exactly their addition `A + B`
     *  (union of nonzero values), which for the Range group is computed by
     *  `range_spadd`. Since `range_spadd` requires operands and result to be in
     *  (device) memory, this overload is selected only when neither input nor
     *  output is buffered (out-of-core).
     */
    template< typename TMutableCRSMatrix, typename TCRSMatrix >
    inline TMutableCRSMatrix
    merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                          crs_matrix::RangeGroup /* group tag */,
                          std::false_type /* is_buffered */ )
    {
      using execution_space = Kokkos::DefaultExecutionSpace;
      execution_space space;

      auto a_entries = dindex1.entries_device_view( space );
      auto a_rowmap = dindex1.rowmap_device_view( space );
      auto b_entries = dindex2.entries_device_view( space );
      auto b_rowmap = dindex2.rowmap_device_view( space );
      auto c_entries = TMutableCRSMatrix::make_entries_device_view( space );
      auto c_rowmap = TMutableCRSMatrix::make_rowmap_device_view( space );

      SparseRangeHandle handle( dindex1, dindex2, space );
      range_spadd( handle, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap,
                   c_entries );

      auto c_nnz = crs_matrix::nnz( c_entries, c_rowmap, crs_matrix::RangeGroup{} );
      return TMutableCRSMatrix( dindex1.numCols(), c_entries, c_rowmap, c_nnz );
    }

    /**
     *  @brief  Merge two distance indices (Range Group)
     *
     *  @param  dindex1 first distance index
     *  @param  dindex2 second distance index
     *  @return a mutable merged distance index of type `TMutableCRSMatrix`
     *
     *  Selects the merge strategy by operand residency: `range_spadd` for in-core
     *  operands, or the host-sequential streaming merge when an input or the
     *  output is buffered (out-of-core).
     */
    template< typename TMutableCRSMatrix, typename TCRSMatrix >
    inline TMutableCRSMatrix
    merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                          crs_matrix::RangeGroup /* tag */ )
    {
      using buffered = std::integral_constant< bool,
          crs_matrix::is_buffered< TCRSMatrix >::value ||
          crs_matrix::is_buffered< TMutableCRSMatrix >::value >;
      return merge_distance_index< TMutableCRSMatrix >(
          dindex1, dindex2, crs_matrix::RangeGroup{}, buffered{} );
    }

    /**
     *  @brief  Merge two distance indices (Basic Group)
     *
     *  @param  dindex1 first distance index
     *  @param  dindex2 second distance index
     *  @return a mutable merged distance index of type `TMutableCRSMatrix`
     *
     *  NOTE: The resulting mutable matrix can be assigned to a immutable compressed
     *        matrix afterwards.
     *
     *  NOTE: The input distance indices are passed by non-const references, since
     *        containers in const Buffered specialisations cannot be iterated.
     */
    template< typename TMutableCRSMatrix, typename TCRSMatrix >
    inline TMutableCRSMatrix
    merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2,
                          crs_matrix::BasicGroup /* tag */ )
    {
      typedef TMutableCRSMatrix crsmat_mutable_type;
      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;
      typedef typename crsmat_type::size_type size_type;

      auto entries = crsmat_mutable_type::make_entries();
      auto rowmap = crsmat_mutable_type::make_rowmap();
      size_type cursor1 = 0;    // current entry index in the first distance index
      size_type cursor2 = 0;    // current entry index in the second distance index
      size_type end1;  // last entry index of the row in the first distance index
      size_type end2;  // last entry index of the row in the second distance index

      assert( dindex1.numRows() == dindex2.numRows() );
      assert( dindex1.numCols() == dindex2.numCols() );

      auto nof_rows = dindex1.numRows();
      auto nof_cols = dindex1.numCols();

      for ( ordinal_type nrow = 0; nrow < nof_rows; ++nrow ) {
        rowmap.push_back( entries.size() );

        end1 = dindex1.rowMap( nrow + 1 );
        end2 = dindex2.rowMap( nrow + 1 );
        while ( cursor1 < end1 ) {
          if ( cursor2 >= end2 ) {
            while ( cursor1 < end1 ) entries.push_back( dindex1.entry( cursor1++ ) );
            break;
          }
          if ( dindex2.entry( cursor2 ) < dindex1.entry( cursor1 ) ) {
            entries.push_back( dindex2.entry( cursor2++ ) );
          } else {
            entries.push_back( dindex1.entry( cursor1 ) );
            if ( dindex2.entry( cursor2 ) == dindex1.entry( cursor1 ) ) ++cursor2;
            ++cursor1;
          }
        }
        while ( cursor2 < end2 ) entries.push_back( dindex2.entry( cursor2++ ) );
      }
      rowmap.push_back( entries.size() );
      assert( cursor1 == dindex1.nnz() );
      assert( cursor2 == dindex2.nnz() );

      return crsmat_mutable_type( nof_cols, std::move( entries ), std::move( rowmap ) );
    }

    template< typename TMutableCRSMatrix, typename TCRSMatrix >
    inline TMutableCRSMatrix
    merge_distance_index( TCRSMatrix& dindex1, TCRSMatrix& dindex2 )
    {
      static_assert(
        crs_matrix::is_same_group< typename TMutableCRSMatrix::spec_type,
                                   typename TCRSMatrix::spec_type >::value,
        "input and output distance indices must be in the same group" );
      return merge_distance_index< TMutableCRSMatrix >(
        dindex1, dindex2, typename crs_matrix::Group< typename TCRSMatrix::spec_type >::type{} );
    }

    /**
     *  @brief  Compress a distance index by removing intra-node loci pairs (Basic Group)
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
    compress_distance_index( TCRSMatrix& dindex, TGraph const& graph,
                             crs_matrix::BasicGroup /* tag */ )
    {
      typedef TGraph graph_type;
      typedef typename graph_type::rank_type rank_type;

      typedef TMutableCRSMatrix crsmat_mutable_type;
      typedef TCRSMatrix crsmat_type;
      typedef typename crsmat_type::ordinal_type ordinal_type;
      typedef typename crsmat_type::size_type size_type;

      auto entries = crsmat_mutable_type::make_entries();
      auto rowmap = crsmat_mutable_type::make_rowmap();
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

    /**
     *  @brief  Compress a distance index by removing intra-node loci pairs (Range Group)
     *
     *  NOTE: Not implemented. Actually, it might not even make sense to
     *  implement this for the Range group, since the compression gain is likely
     *  to be negligible.
     */
    template< typename TMutableCRSMatrix,
              typename TCRSMatrix,
              typename TGraph >
    inline TMutableCRSMatrix
    compress_distance_index( TCRSMatrix& /* dindex */, TGraph const& /* graph */,
                             crs_matrix::RangeGroup /* tag */ )
    {
      static_assert( crs_matrix::always_false< TCRSMatrix >::value,
                     "compress_distance_index is not implemented for the Range group" );
    }

    /**
     *  @brief  Compress a distance index by removing intra-node loci pairs
     *
     *  @param  dindex input distance index
     *  @param  graph underlying graph
     *  @return a mutable compressed distance index of type `TMutableCRSMatrix`
     *
     *  Dispatches to the group-specific implementation inferred from the input
     *  matrix. The output mutable matrix must be in the same group as the input.
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
      static_assert(
        crs_matrix::is_same_group< typename TMutableCRSMatrix::spec_type,
                                   typename TCRSMatrix::spec_type >::value,
        "input and output distance indices must be in the same group" );
      return compress_distance_index< TMutableCRSMatrix >(
        dindex, graph,
        typename crs_matrix::Group< typename TCRSMatrix::spec_type >::type{} );
    }

    /**
     *  @brief  Create a distance index for a graph by its adjacency matrix.
     *
     *  @param  ra     The adjacency matrix of the graph in CRS Range foramt.
     *  @param  dlo    The lower bound of distance constraint.
     *  @param  dup    The upper bound of distance constraint.
     *  @return A RCRS matrix representing the distance index.
     *
     *  NOTE: This method assumes that the input graph is sorted such that node
     *  rank ranges in components are disjoint.
     *
     *  NOTE: The input adjacency matrix should be in Range CRS format.
     */
    template< typename TRCRSMatrix,
              typename TSparseConfig = DefaultSparseConfiguration,
              typename TTimer = Kokkos::Timer >
    inline TRCRSMatrix
    create_distance_index( TRCRSMatrix const& ra, unsigned int dlo,
                           unsigned int dup, TSparseConfig config = {},
                           /* measures whole body minus host<->device copies */
                           TTimer* timer1_ptr = nullptr,
                           /* measures individual operations */
                           TTimer* timer2_ptr = nullptr )
    {
      using rcrsmatrix_t = TRCRSMatrix;
      using config_type = TSparseConfig;
      using execution_space = typename config_type::execution_space;
      using handle_t = SparseRangeHandle< TRCRSMatrix, execution_space >;
      using size_type = typename rcrsmatrix_t::size_type;

      static_assert( crs_matrix::is_range_crs_matrix< rcrsmatrix_t >::value,
                     "matrix should be in Range CRS format." );

      DIVERG_ASSERT( dlo <= dup && dup != 0 );
      assert( ra.numCols() == ra.numRows() );

      execution_space space;

      auto nrows = ra.numRows();
      auto ncols = ra.numCols();

      // rA^dlo . (rA + rI)^(dup - dlo) (device)
      auto d_entries = rcrsmatrix_t::make_entries_device_view( space );
      auto d_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );
      size_type d_nnz = 0;

      {
        // (rA + rI)^(dup - dlo) (device)
        auto raid_entries = rcrsmatrix_t::make_entries_device_view( space );
        auto raid_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );
        size_type raid_nnz = 0;

        // rA^dlo (device)
        auto rad_entries = rcrsmatrix_t::make_entries_device_view( space );
        auto rad_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );
        size_type rad_nnz = 0;

        {
          // rA + rI (device)
          auto rai_entries = rcrsmatrix_t::make_entries_device_view( space );
          auto rai_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );
          size_type rai_nnz = 0;

          {
#ifdef DIVERG_STATS
            if ( timer2_ptr ) timer2_ptr->reset();
#endif
            // rA (device)
            auto ra_entries = ra.entries_device_view( space );
            auto ra_rowmap = ra.rowmap_device_view( space );

#ifdef DIVERG_STATS
            if ( timer2_ptr ) {
                space.fence();
                auto duration = timer2_ptr->seconds();
                std::cout << "diverg::create_distance_index::copy_to_device"
                             " time: " << duration * 1000 << "ms" << std::endl;
            }
#endif

            {
              // rI (device)
              auto i_entries = rcrsmatrix_t::make_entries_device_view( space );
              auto i_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );
#ifdef DIVERG_STATS
              if ( timer1_ptr ) timer1_ptr->reset();
              if ( timer2_ptr ) timer2_ptr->reset();
#endif
              // Computing rI
              create_range_identity_matrix( i_rowmap, i_entries, nrows );
#ifdef DIVERG_STATS
              if ( timer2_ptr ) {
                space.fence();
                auto duration = timer2_ptr->seconds();
                std::cout << "diverg::create_distance_index::create_range_identity_matrix"
                             " time: " << duration * 1000 << "ms" << std::endl;
              }
#endif
              handle_t handle( ncols, nrows, ncols, ra.nnz() ); // rI + rA
#ifdef DIVERG_STATS
              if ( timer2_ptr ) timer2_ptr->reset();
#endif
              // Computing rI + rA
              range_spadd( handle, i_rowmap, i_entries, ra_rowmap, ra_entries,
                           rai_rowmap, rai_entries );
              // rai_nnz = crs_matrix::nnz( rai_entries, rai_rowmap,
              //                            crs_matrix::RangeGroup{} );
              rai_nnz = ra.nnz() + nrows;  // estimate
#ifdef DIVERG_STATS
              if ( timer2_ptr ) {
                space.fence();
                auto duration = timer2_ptr->seconds();
                std::cout << "diverg::create_distance_index::range_spadd time: "
                          << duration * 1000 << "ms" << std::endl;
              }
#endif
            }  // free: rI

            if ( dlo != 0 ) {
              handle_t handle( ncols, ra.nnz(), ncols, ra.nnz() );  // rA^dlo
              // Computing rA^dlo
              rad_nnz = range_power_inplace( handle, ra_rowmap, ra_entries,
                                             rad_rowmap, rad_entries, dlo,
                                             config, timer2_ptr );
            }
          }  // free: rA

          handle_t handle( ncols, rai_nnz, ncols, rai_nnz );  // (rA + rI)^(dup - dlo)
          // Computing (rA + rI)^(dup - dlo)
          raid_nnz = range_power_inplace( handle, rai_rowmap, rai_entries,
                                          raid_rowmap, raid_entries, dup - dlo,
                                          config, timer2_ptr );
        }  // free: (rA + rI)

        if ( dlo != 0 ) {
          handle_t handle( ncols, rad_nnz, ncols, raid_nnz );  // d
          // Computing d
          d_nnz = range_spgemm( handle, rad_rowmap, rad_entries, raid_rowmap,
                                raid_entries, d_rowmap, d_entries, config );
        }
        else {
          d_entries = raid_entries;
          d_rowmap = raid_rowmap;
          d_nnz = raid_nnz;
        }

#ifdef DIVERG_STATS
        if ( timer1_ptr ) {
          space.fence();
          auto duration = timer1_ptr->seconds();
          std::cout << "diverg::create_distance_index time: " << duration * 1000
                    << "ms" << std::endl;
        }
#endif
      } // free: rA^dlo and (rA + rI)^(dup - dlo)
      return TRCRSMatrix( ncols, d_entries, d_rowmap, d_nnz );
    }
  } // namespace util
} // namespace diverg

#endif // DIVERG_DINDEX_HPP__
