/**
 *    @file  range_sparse_utils.hpp
 *   @brief  Utilities for range sparse matrix data structures.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Sat May 04, 2024  12:48
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_RANGE_SPARSE_UTILS_HPP__
#define DIVERG_RANGE_SPARSE_UTILS_HPP__

#include <iostream>
#include <iomanip>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <KokkosKernels_Utils.hpp>

#include "crs_matrix.hpp"
#include "range_sparse_base.hpp"
#include "basic_types.hpp"
#include "utils.hpp"
#include "random.hpp"


namespace diverg {
  /**
   *  @brief  Print a matrix of type `KokkosSparse::CrsMatrix`.
   *
   *  NOTE: All Views associated with the input matrix are assumed to be on
   *  device memory space. If they are host accessible deep copy does not copy
   *  anything.
   */
  template< typename TXCRSMatrix >
  static void
  print( const TXCRSMatrix& m, bool verbose=true, bool print_all=false )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type::HostMirror const_host_entries_t;
    typedef typename xcrsmatrix_t::values_type::non_const_type::HostMirror const_host_values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type::HostMirror const_host_row_map_t;
    typedef typename xcrsmatrix_t::ordinal_type ordinal_type;

    auto label = m.values.label();
    auto nrows = m.numRows();
    auto ncols = m.numCols();
    auto nnz = m.nnz();

    std::cout << "[INFO] Matrix '" << label << "'"
              << " (" << nrows << "x" << ncols << ") with " << nnz << " non-zero elements:\n";

    if ( verbose ) {
      const_host_entries_t entries = Kokkos::create_mirror_view( m.graph.entries );
      const_host_values_t values = Kokkos::create_mirror_view( m.values );
      const_host_row_map_t rowmap = Kokkos::create_mirror_view( m.graph.row_map );

      Kokkos::deep_copy( entries, m.graph.entries );
      Kokkos::deep_copy( values, m.values );
      Kokkos::deep_copy( rowmap, m.graph.row_map );

      std::cout << "   ... ┬─\n";
      std::cout << "   ... ├─ entries  (" << entries.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( entries, print_all );
      std::cout << "   ... ├─ values   (" << values.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( values, print_all );
      std::cout << "   ... ╰─ row map  (" << rowmap.extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( rowmap, print_all );
      std::cout << "   ... " << std::endl;

      auto width = std::log( nrows ) + 1;
      std::cout << "   ... " << std::setw( width - 2 ) << label << " = [" << std::endl;
      for ( ordinal_type i = 0; i < nrows; ++i ) {
        std::cout << "   ... " << std::setw( width ) << i << ": ";
        auto end = rowmap( i + 1 );
        for ( auto j = rowmap( i ); j < end; ++j ) {
          std::cout << " " << entries( j );
        }
        std::cout << "\n";
      }
      std::cout << "   ... " << std::setw( width + 2 ) << "]" << std::endl;
    }
  }

  /**
   *  @brief  Print a matrix of type 'Range' spacialised `diverg::CRSMatrix`.
   *
   *  NOTE: All Views associated with the matrix are assumed to be on host
   *  memory space.
   */
  template< typename TSpec, typename TOrdinal, typename TSize >
  static /* void */ std::enable_if_t<
      std::is_same< typename diverg::crs_matrix::Group< TSpec >::type,
                    diverg::crs_matrix::RangeGroup >::value >
  print( diverg::CRSMatrix< TSpec, bool, TOrdinal, TSize >& m,
         std::string label = {}, bool verbose = true, bool print_all = false )
  {
    typedef TOrdinal ordinal_type;

    if ( label.empty() ) label = "A";

    std::cout << "[INFO] Matrix '" << label << "'"
              << " (" << m.numRows() << "x" << m.numCols() << ") with "
              << m.nnz() << " non-zero elements:\n";

    if ( verbose ) {
      std::cout << "   ... ┬─\n";
      std::cout << "   ... ├─ entries  (" << m.entries_view().extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( m.entries_view(), print_all );
      std::cout << "   ... ├─ values   (" << m.entries_view().extent( 0 ) << "): 1 ... 1\n";
      std::cout << "   ... ╰─ row map  (" << m.rowmap_view().extent( 0 ) << "): ";
      KokkosKernels::Impl::print_1Dview( m.rowmap_view(), print_all );
      std::cout << "   ... " << std::endl;

      auto width = std::log( m.numRows() ) + 1;
      std::cout << "   ... " << std::setw( width - 2 ) << label << " = [" << std::endl;
      for ( ordinal_type i = 0; i < m.numRows(); ++i ) {
        std::cout << "   ... " << std::setw( width ) << i << ": ";
        auto end = m.rowMap( i + 1 );
        for ( auto j = m.rowMap( i ); j < end; j += 2 ) {
          std::cout << " (" << m.entry( j ) << ", " << m.entry( j + 1 ) << ")";
        }
        std::cout << "\n";
      }
      std::cout << "   ... " << std::setw( width + 2 ) << "]" << std::endl;
    }
  }

  /**
   *  @brief  Create a random square matrix (on-host) [SLOW].
   *
   *  @param  n     order of the output matrix
   *  @param  nnz   number of non-zero values
   *  @param  lower lower bound of value range (inclusive)
   *  @param  upper upper bound of value range (exclusive)
   *
   *  @return a `KokkosSparse::CrsMatrix`-like square matrix of order `n` with
   *  `nnz` number of non-zero random values in range [lower, upper).
   *
   *  NOTE: Use `create_random_matrix` instead.
   *
   *  NOTE: The output matrix is constructed on **host** memory space, then it
   *  is transferred (deep-copied) to device memory space. It is suitable for
   *  small matrices.
   *
   *  NOTE: The final matrix is on device memory and host mirrors get
   *  deallocated on return.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_random_matrix_on_host( typename TXCRSMatrix::ordinal_type n,
                                typename TXCRSMatrix::size_type nnz,
                                typename TXCRSMatrix::value_type lower=std::numeric_limits< typename TXCRSMatrix::value_type >::min(),
                                typename TXCRSMatrix::value_type upper=std::numeric_limits< typename TXCRSMatrix::value_type >::max() )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::value_type value_type;
    typedef typename xcrsmatrix_t::size_type size_type;
    typedef typename xcrsmatrix_t::values_type::non_const_type values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type row_map_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type entries_t;

    assert( n > 1 && nnz > 0 && ( nnz / n ) <= static_cast< size_type >( n ) );

    values_t a_values( Kokkos::ViewAllocateWithoutInitializing( "R" ), nnz );
    row_map_t a_row_map( Kokkos::ViewAllocateWithoutInitializing( "rowmap" ),
                         n + 1 );
    entries_t a_entries( Kokkos::ViewAllocateWithoutInitializing( "entries" ),
                         nnz );

    auto h_a_entries = Kokkos::create_mirror_view( a_entries );
    auto h_a_values = Kokkos::create_mirror_view( a_values );
    auto h_a_row_map = Kokkos::create_mirror_view( a_row_map );

    // Zero initialisation: the rest will be initialized later
    h_a_row_map( 0 ) = 0;

    Kokkos::parallel_for(
        "diverg::crs_matrix::::create_random_matrix_on_host::random_values",
        Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, nnz ),
        [=]( const uint64_t i ) {
          value_type v = 0;
          while ( v == 0 ) v = diverg::random::random_integer( lower, upper );
          h_a_values( i ) = v;
        } );

    {
      // Distributing nnz values into rows
      std::size_t i = 0;
      while ( i < nnz ) {
        auto idx = diverg::random::random_index( n );
        do {
          if ( h_a_row_map( idx + 1 ) < static_cast< size_type >( n ) ) {
            ++h_a_row_map( idx + 1 );
            ++i;
            break;
          }
          idx = ( idx + 1 ) % n;
        } while( true );
      }

      //for ( i = 1; i < n; ++i ) h_a_row_map( i + 1 ) += h_a_row_map( i );
      Kokkos::parallel_scan(
          "diverg::crs_matrix::::create_random_matrix_on_host::compute_row_map",
          Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, n ),
          [=]( const int i, size_type& update, const bool final ) {
            // Load old value in case we update it before accumulating
            const size_type val_ip1 = h_a_row_map( i + 1 );
            update += val_ip1;
            if ( final )
              h_a_row_map( i + 1 )
                  = update;  // only update array on final pass
          } );
    }

    Kokkos::parallel_for(
        "diverg::crs_matrix::::create_random_matrix_on_host::random_entries",
        Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >( 0, n ),
        [=]( const uint64_t i ) {
          auto l = h_a_row_map( i );
          auto u = h_a_row_map( i + 1 );
          auto begin = h_a_entries.data() + l;
          auto end = h_a_entries.data() + u;
          std::sample( diverg::RangeIterator< decltype( n ) >{ 0 },
                       diverg::RangeIterator{ n }, begin, u - l,
                       diverg::random::gen );
          std::sort( begin, end );
        } );

    Kokkos::deep_copy( a_entries, h_a_entries );
    Kokkos::deep_copy( a_values, h_a_values );
    Kokkos::deep_copy( a_row_map, h_a_row_map );

    assert( h_a_row_map( n ) == nnz );

    return xcrsmatrix_t( "Random Matrix", n, n, nnz, a_values, a_row_map, a_entries );
  }

  /**
   *  @brief  Create a random square matrix (on-device) [FAST].
   *
   *  @param  n     order of the output matrix
   *  @param  nnz   number of non-zero values
   *  @param  lower lower bound of value range (inclusive)
   *  @param  upper upper bound of value range (exclusive)
   *
   *  @return a `KokkosSparse::CrsMatrix`-like square matrix of order `n` with
   *  `nnz` number of non-zero random values in range [lower, upper).
   *
   *  This function is usually much faster than `create_random_matrix_on_host`.
   *
   *  NOTE: The output matrix is constructed on device memory space. It does
   *  not performs any deep copy between device and host.
   */
  template< typename TXCRSMatrix,
            typename THostSpace=Kokkos::DefaultHostExecutionSpace >
  inline TXCRSMatrix
  create_random_matrix(
      typename TXCRSMatrix::ordinal_type n,
      typename TXCRSMatrix::size_type nnz,
      typename TXCRSMatrix::value_type lower
      = std::numeric_limits< typename TXCRSMatrix::value_type >::min(),
      typename TXCRSMatrix::value_type upper
      = std::numeric_limits< typename TXCRSMatrix::value_type >::max() )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename TXCRSMatrix::execution_space execution_space;
    typedef typename xcrsmatrix_t::value_type value_type;
    typedef typename xcrsmatrix_t::ordinal_type ordinal_type;
    typedef typename xcrsmatrix_t::size_type size_type;
    typedef typename xcrsmatrix_t::values_type::non_const_type values_t;
    typedef typename xcrsmatrix_t::row_map_type::non_const_type row_map_t;
    typedef typename xcrsmatrix_t::index_type::non_const_type entries_t;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;
    typedef typename Kokkos::Random_XorShift64_Pool< execution_space > random_pool_t;
    typedef typename random_pool_t::generator_type generator_t;

    assert( n > 1 && nnz > 0 && ( nnz / n ) <= static_cast< size_type >( n ) );

    // NOTE: (nnz/n) should fit in `ordinal_type` (asserted above) as
    //       `sizeof( size_type ) >= sizeof( ordinal_type )`
    // NOTE: `ordinal_type` to match with the signness of `n`
    ordinal_type nnz_per_row = nnz / n;

    values_t r_values( Kokkos::ViewAllocateWithoutInitializing( "R" ), nnz );
    row_map_t r_row_map( "rowmap", n + 1 );
    entries_t r_entries( Kokkos::ViewAllocateWithoutInitializing( "entries" ),
                         nnz );

    random_pool_t random_pool( diverg::random::rd() );

    Kokkos::parallel_for( "diverg::crs_matrix::create_random_matrix::random_values",
                          range_policy_t( 0, nnz ),
                          KOKKOS_LAMBDA ( const uint64_t i ) {
                            generator_t generator = random_pool.get_state();
                            value_type v = 0;
                            while ( v == 0 ) {
                              v = Kokkos::rand< generator_t,
                                                value_type >::draw( generator,
                                                                    lower,
                                                                    upper );
                            }
                            r_values( i ) = v;
                            random_pool.free_state( generator );
                          } );

    if ( nnz_per_row == n ) {  // if nnz = n*n (n*n may be really large to fit in an integer)
      Kokkos::parallel_for(
          "diverg::crs_matrix::create_random_matrix::fill_nnz",
          range_policy_t( 0, n ),
          KOKKOS_LAMBDA ( const uint64_t i ) { r_row_map( i + 1 ) = n; } );
    }
    else {
      auto d_nnz = nnz;
      if ( nnz_per_row > n / 2 ) {  // in this case, distribute the zero values is cheaper.
        d_nnz = ( n - ( nnz_per_row ) ) * n - ( nnz % n );
      }

      Kokkos::parallel_for(
          "diverg::crs_matrix::create_random_matrix::distribute_nnz",
          range_policy_t( 0, d_nnz ), KOKKOS_LAMBDA ( const uint64_t ) {
            generator_t generator = random_pool.get_state();

            ordinal_type idx =
              Kokkos::rand< generator_t, ordinal_type >::draw( generator, n );

            bool exchanged = false;
            do {
              auto ptr = &r_row_map( idx + 1 );
              auto value = Kokkos::atomic_load( ptr );
              while ( value < static_cast< size_type >( n ) ) {
                exchanged = ( Kokkos::atomic_compare_exchange(
                    ptr, value, value + 1 ) == value );
                if ( exchanged ) break;
                value = Kokkos::atomic_load( ptr );
              }
              idx = ( idx + 1 ) % n;
            } while ( !exchanged );

            random_pool.free_state( generator );
          } );

      if ( d_nnz != nnz ) {
        Kokkos::parallel_for(
            "diverg::crs_matrix::create_random_matrix::reverse_nnz_dist",
            range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t i ) {
              r_row_map( i + 1 ) = n - r_row_map( i + 1 );
            } );
      }
    }

    Kokkos::parallel_scan(
        "diverg::crs_matrix::create_random_matrix::compute_row_map",
        range_policy_t( 0, n ),
        KOKKOS_LAMBDA ( const int i, size_type& partial_sum,
                        const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type value = r_row_map( i + 1 );
          partial_sum += value;
          // only update array on final pass
          if ( final ) r_row_map( i + 1 ) = partial_sum;
          if ( i == 0 ) r_row_map( 0 ) = 0;
        } );

    Kokkos::parallel_for(
        "diverg::crs_matrix::create_random_matrix::random_entries",
        range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t i ) {
          auto l = r_row_map( i );
          auto u = r_row_map( i + 1 );
          auto begin = r_entries.data() + l;
          assert( l <= u );
          ordinal_type k = u - l;
          if ( k != 0 ) {
            generator_t generator = random_pool.get_state();

            // Reservoir sampling algorithm
            ordinal_type j = 0;
            for ( j = 0; j < k; ++j ) *( begin + j ) = j;
            for ( ; j < n; ++j ) {
              ordinal_type r = Kokkos::rand< generator_t, ordinal_type >::draw(
                  generator, j );
              if ( r < k ) *( begin + r ) = j;
            }

            random_pool.free_state( generator );
          }
        } );

    auto func = SortEntriesFunctor( r_row_map, r_entries );

    Kokkos::parallel_for(
        "diverg::crs_matrix::create_random_matrix::sort_entries",
        func.policy( n ), func );

    return xcrsmatrix_t( "Random Matrix", n, n, nnz, r_values, r_row_map,
                         r_entries );
  }

  /**
   *  @brief  Create a random binary `KokkosSparse::CrsMatrix`.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_random_binary_matrix( typename TXCRSMatrix::ordinal_type n,
                               typename TXCRSMatrix::size_type nnz )
  {
    return create_random_matrix< TXCRSMatrix >( n, nnz, 1, 2 );
  }

  /**
   *  @brief  Create a random binary matrix and return in both KokkosKernels
   *          CRS and DIVERG Range CRS formats.
   */
  template< typename TXCRSMatrix, typename TRCRSMatrix >
  inline TXCRSMatrix
  create_random_binary_matrix( typename TXCRSMatrix::ordinal_type n,
                               typename TXCRSMatrix::size_type nnz,
                               TRCRSMatrix& range_crs )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename TXCRSMatrix::HostMirror xcrsmatrix_host_t;
    typedef TRCRSMatrix range_crsmatrix_t;

    typename xcrsmatrix_t::index_type::non_const_type::HostMirror h_r_entries;
    typename xcrsmatrix_t::values_type::non_const_type::HostMirror h_r_values;
    typename xcrsmatrix_t::row_map_type::non_const_type::HostMirror h_r_row_map;

    auto crs = create_random_binary_matrix< xcrsmatrix_t >( n, nnz );

    h_r_values = Kokkos::create_mirror_view( crs.values );
    h_r_row_map = Kokkos::create_mirror_view( crs.graph.row_map );
    h_r_entries = Kokkos::create_mirror_view( crs.graph.entries );

    Kokkos::deep_copy( h_r_values, crs.values );
    Kokkos::deep_copy( h_r_row_map, crs.graph.row_map );
    Kokkos::deep_copy( h_r_entries, crs.graph.entries );

    const xcrsmatrix_host_t h_crs( "R host copy", n, n, nnz, h_r_values, h_r_row_map, h_r_entries );
    range_crs = range_crsmatrix_t( h_crs );

    return crs;
  }
}  /* -----  end of namespace diverg  ----- */


#endif // DIVERG_RANGE_SPARSE_UTILS_HPP__
