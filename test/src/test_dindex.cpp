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

  auto a_entries = actual.entries_view();
  auto e_entries = expected.entries_view();
  REQUIRE( a_entries.extent( 0 ) == e_entries.extent( 0 ) );
  for ( std::size_t i = 0; i < e_entries.extent( 0 ); ++i )
    REQUIRE( a_entries( i ) == e_entries( i ) );

  auto a_rowmap = actual.rowmap_view();
  auto e_rowmap = expected.rowmap_view();
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
    std::string graph_path = test_data_dir + "/middle/m.gfa";
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
  }
}
