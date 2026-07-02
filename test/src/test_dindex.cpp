/**
 *    @file  test_dindex.cpp
 *   @brief  Test scenarios for the distance index module (`dindex.hpp`).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Sun Jun 21, 2026
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2026, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdint>
#include <string>

#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <diverg/dindex.hpp>
#include <diverg/dindex_eti.hpp>

#include <KokkosSparse_CrsMatrix.hpp>

#include "test_base.hpp"


using namespace diverg;

/**
 *  @brief  Assert that two Range CRS matrices are byte-for-byte identical.
 */
template< typename TRangeCRSMatrix >
inline void
check_identical( TRangeCRSMatrix& actual, TRangeCRSMatrix& expected )
{
  REQUIRE( actual.numRows() == expected.numRows() );
  REQUIRE( actual.numCols() == expected.numCols() );
  REQUIRE( actual.nnz() == expected.nnz() );

  auto a_entries = diverg::entries_view( actual );
  auto e_entries = diverg::entries_view( expected );
  REQUIRE( a_entries.extent( 0 ) == e_entries.extent( 0 ) );
  for ( std::size_t i = 0; i < e_entries.extent( 0 ); ++i )
    REQUIRE( a_entries( i ) == e_entries( i ) );

  auto a_rowmap = diverg::rowmap_view( actual );
  auto e_rowmap = diverg::rowmap_view( expected );
  REQUIRE( a_rowmap.extent( 0 ) == e_rowmap.extent( 0 ) );
  for ( std::size_t i = 0; i < e_rowmap.extent( 0 ); ++i )
    REQUIRE( a_rowmap( i ) == e_rowmap( i ) );
}

SCENARIO( "Graph adjacency matrix construction", "[dindex]" )
{
  typedef int scalar_t;
  typedef Kokkos::DefaultHostExecutionSpace host_space;
  typedef KokkosSparse::CrsMatrix< scalar_t, int32_t, host_space > xcrsmatrix_t;
  typedef diverg::CRSMatrix< diverg::crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > range_crsmatrix_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A sequence graph" )
  {
    std::string graph_path = test_data_dir + "/multi/multi.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );

    WHEN( "Range CRS adjacency matrix is directly constructed for whole-graph" )
    {
      auto actual = diverg::util::range_adjacency_matrix< range_crsmatrix_t >( graph );

      THEN( "It should be identical to the adjacency matrix in rCRS converted from CRS" )
      {
        auto a = diverg::util::adjacency_matrix< xcrsmatrix_t >( graph );
        range_crsmatrix_t expected( a );
        REQUIRE( expected.numRows() > 0 );
        REQUIRE( expected.nnz() > 0 );
        check_identical( actual, expected );
      }
    }

    WHEN( "Range CRS adjacency matrix is directly constructed for components" )
    {
      THEN( "It should be identical to the adjacency matrix in rCRS converted from CRS" )
      {
        typename graph_type::rank_type lower = 1;
        auto node_count = graph.get_node_count();
        std::size_t nof_checked = 0;

        while ( lower <= node_count ) {
          auto upper = std::get< 0 >(
              diverg::util::upper_node_rank_in_component( graph, lower ) );

          auto a = diverg::util::adjacency_matrix< xcrsmatrix_t >(
              graph, lower, upper );
          range_crsmatrix_t expected( a );
          auto actual = diverg::util::range_adjacency_matrix< range_crsmatrix_t >(
              graph, lower, upper );

          check_identical( actual, expected );

          ++nof_checked;
          lower = upper;
        }

        REQUIRE( nof_checked > 0 );
      }
    }

    WHEN( "A block-diagonal range matrix is assembled from per-component Range blocks" )
    {
      typedef range_crsmatrix_t::ordinal_type ordinal_t;

      // Reference: the whole-graph range adjacency.
      auto expected = diverg::util::range_adjacency_matrix< range_crsmatrix_t >( graph );

      // Assemble the same matrix from per-component Range blocks.
      auto nrows = static_cast< ordinal_t >( gum::util::total_nof_loci( graph ) );
      auto provider = [&graph]( auto partial ) {
        typename graph_type::rank_type lower = 1;
        auto node_count = graph.get_node_count();
        while ( lower <= node_count ) {
          auto upper = std::get< 0 >(
              diverg::util::upper_node_rank_in_component( graph, lower ) );
          auto block
              = diverg::util::range_adjacency_matrix< range_crsmatrix_t >(
                  graph, lower, upper );
          auto soff = static_cast< ordinal_t >(
              gum::util::id_to_charorder( graph, graph.rank_to_id( lower ) ) );
          partial( block, soff, soff );
          lower = upper;
        }
      };
      range_crsmatrix_t actual( nrows, nrows, provider );

      THEN( "It should be identical to the directly-built whole-graph range matrix" )
      {
        check_identical( actual, expected );
      }
    }

    WHEN( "A block-diagonal range matrix is assembled from per-component Basic blocks" )
    {
      typedef diverg::CRSMatrix< diverg::crs_matrix::Dynamic, bool, uint32_t, uint64_t > basic_crsmatrix_t;
      typedef range_crsmatrix_t::ordinal_type ordinal_t;

      // Reference: the whole-graph range adjacency.
      auto expected = diverg::util::range_adjacency_matrix< range_crsmatrix_t >( graph );

      // Assemble from per-component *native* Basic blocks.
      auto nrows = static_cast< ordinal_t >( gum::util::total_nof_loci( graph ) );
      auto provider = [&graph]( auto partial ) {
        typename graph_type::rank_type lower = 1;
        auto node_count = graph.get_node_count();
        while ( lower <= node_count ) {
          auto upper = std::get< 0 >(
              diverg::util::upper_node_rank_in_component( graph, lower ) );
          auto a = diverg::util::adjacency_matrix< xcrsmatrix_t >( graph, lower, upper );
          basic_crsmatrix_t block( a );  // KokkosSparse -> diverg Basic
          auto soff = static_cast< ordinal_t >(
              gum::util::id_to_charorder( graph, graph.rank_to_id( lower ) ) );
          partial( block, soff, soff );
          lower = upper;
        }
      };
      range_crsmatrix_t actual( nrows, nrows, provider );

      THEN( "It should be identical to the directly-built whole-graph range matrix" )
      {
        check_identical( actual, expected );
      }
    }
  }
}

SCENARIO( "Distance index construction templated/compiled boundary", "[dindex]" )
{
  typedef diverg::CRSMatrix< diverg::crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > range_crsmatrix_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A range adjacency matrix and a distance range" )
  {
    std::string graph_path = test_data_dir + "/multi/multi.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );

    unsigned int dmin = 1;
    unsigned int dmax = 8;

    auto ra = diverg::util::range_adjacency_matrix< range_crsmatrix_t >( graph );

    WHEN( "the compiled (non-template) boundary is called" )
    {
      auto boundary = diverg::util::create_distance_index( ra, dmin, dmax );
      auto templated = diverg::util::create_distance_index(
          ra, dmin, dmax, diverg::DefaultSparseConfiguration{} );

      THEN( "it is identical to the templated implementation" )
      {
        REQUIRE( boundary.numRows() > 0 );
        check_identical( boundary, templated );
      }
    }
  }
}

SCENARIO( "create_distance_index ETI supports the uint32/uint32 combination", "[dindex]" )
{
  typedef diverg::CRSMatrix< diverg::crs_matrix::RangeDynamic, bool, uint32_t, uint32_t > range_crsmatrix32_t;
  typedef gum::SeqGraph< gum::Succinct > graph_type;

  GIVEN( "A uint32/uint32 range adjacency matrix" )
  {
    std::string graph_path = test_data_dir + "/multi/multi.gfa";
    graph_type graph;
    gum::util::load( graph, graph_path, gum::util::GFAFormat{}, true );

    auto ra = diverg::util::range_adjacency_matrix< range_crsmatrix32_t >( graph );

    WHEN( "the compiled (uint32/uint32) boundary is called" )
    {
      auto boundary = diverg::util::create_distance_index( ra, 1, 3 );
      auto templated = diverg::util::create_distance_index(
          ra, 1, 3, diverg::DefaultSparseConfiguration{} );

      THEN( "it is identical to the templated implementation" )
      {
        REQUIRE( boundary.numRows() > 0 );
        check_identical( boundary, templated );
      }
    }
  }
}
