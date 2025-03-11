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
      using rcrs_spec_type = typename rcrsmatrix_t::spec_type;
      using group_type = typename crs_matrix::Group< rcrs_spec_type >::type;

      static_assert( std::is_same< group_type, crs_matrix::RangeGroup >::value,
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
