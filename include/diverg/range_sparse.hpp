/**
 *    @file  range_sparse.hpp
 *   @brief  Range sparse matrix operations
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu Sep 07, 2023  20:20
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_RANGE_SPARSE_HPP_
#define DIVERG_RANGE_SPARSE_HPP_

#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <KokkosKernels_Utils.hpp>

#include "crs_matrix.hpp"
#include "range_sparse_base.hpp"
#include "hbitvector.hpp"
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
   * memory space.
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
   *  @brief  Create the identity matrix of order `n` in CRS format.
   *
   *  NOTE: `TXCRSMatrix` should be a `KokkosSparse::CrsMatrix`-like type.
   */
  template< typename TXCRSMatrix >
  inline TXCRSMatrix
  create_identity_matrix( typename TXCRSMatrix::ordinal_type n )
  {
    typedef TXCRSMatrix xcrsmatrix_t;
    typedef typename xcrsmatrix_t::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;

    typename xcrsmatrix_t::values_type::non_const_type i_values(
        Kokkos::ViewAllocateWithoutInitializing( "I" ), n );
    typename xcrsmatrix_t::row_map_type::non_const_type i_row_map(
        Kokkos::ViewAllocateWithoutInitializing( "I_rowmap" ), n + 1 );
    typename xcrsmatrix_t::index_type::non_const_type i_entries(
        Kokkos::ViewAllocateWithoutInitializing( "I_entries" ), n );

    Kokkos::parallel_for(
        "diverg::crs_matrix::create_identity_matrix", range_policy_t( 0, n ),
        KOKKOS_LAMBDA ( const uint64_t i ) {
          i_values( i ) = 1;
          i_row_map( i + 1 ) = i + 1;
          i_entries( i ) = i;
          if ( i == 0 ) i_row_map( 0 ) = 0;
        } );

    return xcrsmatrix_t( "Identity Matrix", n, n, n, i_values, i_row_map, i_entries );
  }

  template< typename TRowMapDeviceView, typename TEntriesDeviceView >
  inline void
  create_range_identity_matrix(
      TRowMapDeviceView& i_rowmap, TEntriesDeviceView& i_entries,
      typename TEntriesDeviceView::value_type n )
  {
    typedef typename TEntriesDeviceView::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > range_policy_t;

    i_entries = TEntriesDeviceView(
        Kokkos::ViewAllocateWithoutInitializing( "I" ), n * 2 );
    i_rowmap = TRowMapDeviceView(
        Kokkos::ViewAllocateWithoutInitializing( "I_rowmap" ), n + 1 );

    Kokkos::parallel_for(
        "diverg::crs_matrix::create_range_identity_matrix",
        range_policy_t( 0, n ), KOKKOS_LAMBDA ( const uint64_t ii ) {
          i_entries( ii * 2 ) = ii;
          i_entries( ii * 2 + 1 ) = ii;
          i_rowmap( ii + 1 ) = ( ii + 1 ) * 2;
          if ( ii == 0 ) i_rowmap( 0 ) = 0;
        } );
  }

  /**
   *  @brief  Create the identity matrix of order `n` in RCRS format.
   *
   *  NOTE: `TRCRSMatrix` should be a `diverg::CRSMatrix`-like type from Range group.
   */
  template< typename TRCRSMatrix, typename TExecSpace=Kokkos::DefaultExecutionSpace >
  inline TRCRSMatrix
  create_range_identity_matrix( typename TRCRSMatrix::ordinal_type n, TExecSpace space={} )
  {
    typedef TRCRSMatrix rcrsmatrix_t;
    typedef typename rcrsmatrix_t::spec_type rcrs_spec_type;
    typedef typename crs_matrix::Group< rcrs_spec_type >::type group_type;

    static_assert( std::is_same< group_type, crs_matrix::RangeGroup >::value,
                   "matrix should be in Range CRS format." );

    auto i_entries = rcrsmatrix_t::make_entries_device_view( space );
    auto i_rowmap = rcrsmatrix_t::make_rowmap_device_view( space );

    create_range_identity_matrix( i_rowmap, i_entries, n );

    return rcrsmatrix_t( n, i_entries, i_rowmap, n /* nnz */ );
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
   * `nnz` number of non-zero random values in range [lower, upper).
   *
   *  This function is usually much faster than `create_random_matrix_on_host`.
   *
   *  NOTE: The output matrix is constructed on device memory space. It does
   *    not performs any deep copy between device and host.
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
                exchanged = Kokkos::atomic_compare_exchange_strong(
                    ptr, value, value + 1 );
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

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC >
  inline void
  range_spadd_symbolic( THandle&,
                        TRowMapDeviceViewA a_rowmap,
                        TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap,
                        TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC& c_rowmap )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spadd_symbolic::count_row_nnz",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          auto b_idx = b_rowmap( row );
          auto b_end = b_rowmap( row + 1 );

          assert( ( a_end - a_idx ) % 2 == 0 );
          assert( ( b_end - b_idx ) % 2 == 0 );

          size_type count = 0;
          if ( a_idx < a_end && b_idx < b_end ) {
            ordinal_type lo;
            ordinal_type hi;
            if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
              lo = a_entries( a_idx );
              hi = a_entries( a_idx + 1 );
              a_idx += 2;
            }
            else {
              lo = b_entries( b_idx );
              hi = b_entries( b_idx + 1 );
              b_idx += 2;
            }

            while ( a_idx < a_end && b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }
              else {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, rs );
                hi = DIVERG_MACRO_MAX( hi, re );
                continue;
              }
              lo = rs;
              hi = re;
              count += 2;
            }

            while ( a_idx < a_end || b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_idx == a_end ) {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }
              else {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, rs );
                hi = DIVERG_MACRO_MAX( hi, re );
                continue;
              }
              count += 2;
              break;
            }
            count += 2;  // last one from the previous loop
          }

          count += ( b_end - b_idx ) + ( a_end - a_idx );  // the rest

          c_rowmap( row + 1 ) = count;
          if ( row == 0 ) c_rowmap( 0 ) = 0;
        } );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spadd_symbolic::computing_row_map_c",
        policy_type( 0, a_nrows ),
        KOKKOS_LAMBDA ( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC >
  inline void
  range_spadd_numeric( THandle&,
                       TRowMapDeviceViewA a_rowmap,
                       TEntriesDeviceViewA a_entries,
                       TRowMapDeviceViewB b_rowmap,
                       TEntriesDeviceViewB b_entries,
                       TRowMapDeviceViewC c_rowmap,
                       TEntriesDeviceViewC& c_entries )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spadd_numeric::count_row_nnz",
        policy_type( 0, a_nrows ), KOKKOS_LAMBDA ( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          auto b_idx = b_rowmap( row );
          auto b_end = b_rowmap( row + 1 );

          assert( ( a_end - a_idx ) % 2 == 0 );
          assert( ( b_end - b_idx ) % 2 == 0 );

          size_type c_idx = c_rowmap( row );
          if ( a_idx < a_end && b_idx < b_end ) {
            ordinal_type lo;
            ordinal_type hi;
            if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
              lo = a_entries( a_idx );
              hi = a_entries( a_idx + 1 );
              a_idx += 2;
            }
            else {
              lo = b_entries( b_idx );
              hi = b_entries( b_idx + 1 );
              b_idx += 2;
            }

            while ( a_idx < a_end && b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_entries( a_idx ) <= b_entries( b_idx ) ) {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }
              else {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, rs );
                hi = DIVERG_MACRO_MAX( hi, re );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = rs;
              hi = re;
            }

            while ( a_idx < a_end || b_idx < b_end ) {
              ordinal_type rs;  // range start
              ordinal_type re;  // range end
              if ( a_idx == a_end ) {
                rs = b_entries( b_idx );
                re = b_entries( b_idx + 1 );
                b_idx += 2;
              }
              else {
                rs = a_entries( a_idx );
                re = a_entries( a_idx + 1 );
                a_idx += 2;
              }

              if ( rs <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, rs );
                hi = DIVERG_MACRO_MAX( hi, re );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = rs;
              hi = re;
              break;
            }
            c_entries( c_idx++ ) = lo;
            c_entries( c_idx++ ) = hi;
          }
          for ( ; a_idx < a_end; ++a_idx )
            c_entries( c_idx++ ) = a_entries( a_idx );
          for ( ; b_idx < b_end; ++b_idx )
            c_entries( c_idx++ ) = b_entries( b_idx );
        } );
  }

  /**
   *  @brief Computing matrix c as the addition of a and b.
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TTimer = Kokkos::Timer >
  inline void
  range_spadd( THandle& handle,
               TRowMapDeviceViewA a_rowmap, TEntriesDeviceViewA a_entries,
               TRowMapDeviceViewB b_rowmap, TEntriesDeviceViewB b_entries,
               TRowMapDeviceViewC& c_rowmap, TEntriesDeviceViewC& c_entries,
               TTimer* timer_ptr = nullptr )
  {
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;

    assert( handle.a_ncols == handle.b_ncols );
    assert( a_rowmap.extent( 0 ) == b_rowmap.extent( 0 ) );

    ordinal_type n = a_rowmap.extent( 0 ) - 1;
    c_rowmap = c_row_map_type(
        Kokkos::ViewAllocateWithoutInitializing( "c_rowmap" ),
        a_rowmap.extent( 0 ) );

#ifdef DIVERG_STATS
    if ( timer_ptr ) timer_ptr->reset();
#endif

    range_spadd_symbolic( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                          c_rowmap );

#ifdef DIVERG_STATS
    if ( timer_ptr ) {
      execution_space{}.fence();
      auto d = timer_ptr->seconds();
      std::cout << "diverg::range_spadd_symbolic time: " << d * 1000 << "ms"
                << std::endl;
    }
#endif

    size_type c_rnnz;
    Kokkos::deep_copy( c_rnnz, Kokkos::subview( c_rowmap, n ) );
    c_entries = c_entries_type( Kokkos::ViewAllocateWithoutInitializing( "C" ),
                                c_rnnz );

#ifdef DIVERG_STATS
    if ( timer_ptr ) timer_ptr->reset();
#endif

    range_spadd_numeric( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                         c_rowmap, c_entries );

#ifdef DIVERG_STATS
    if ( timer_ptr ) {
      execution_space{}.fence();
      auto d = timer_ptr->seconds();
      std::cout << "diverg::range_spadd_numeric time: " << d * 1000 << "ms"
                << std::endl;
    }
#endif
  }

  template< typename TRCRSMatrix,
            typename TExecSpace = Kokkos::DefaultExecutionSpace >
  inline TRCRSMatrix
  range_spadd( TRCRSMatrix const& a, TRCRSMatrix const& b,
               TExecSpace space = {} )
  {
    typedef TRCRSMatrix range_crsmatrix_t;

    assert( a.numCols() == b.numCols() );
    assert( a.numRows() == b.numRows() );

    auto a_entries = a.entries_device_view( space );
    auto a_rowmap = a.rowmap_device_view( space );
    auto b_entries = b.entries_device_view( space );
    auto b_rowmap = b.rowmap_device_view( space );

    auto c_entries = range_crsmatrix_t::make_entries_device_view( space );
    auto c_rowmap = range_crsmatrix_t::make_rowmap_device_view( space );

    SparseRangeHandle handle( a, b, space );

    range_spadd( handle, a_rowmap, a_entries, b_rowmap, b_entries, c_rowmap,
                 c_entries );

    auto nnz = crs_matrix::nnz( c_entries, c_rowmap, crs_matrix::RangeGroup{} );

    // FIXME: since entries and rowmap arrays of range CRS is not a view, there
    // would be an extra copy here and the `c_entries` and `c_rowmap` cannot be
    // moved when the views are on the same memory space (the RCRS ctor does call
    // `deep_copy`).
    return TRCRSMatrix( a.numCols(), c_entries, c_rowmap, nnz );
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (ThreadRangePolicyPartition-BTreeAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of
   *        size `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TExecGrid >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  _range_spgemm_symbolic( THandle&,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          TExecGrid, ThreadRangePolicyPartition, BTreeAccumulator )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;
    typedef std::map< ordinal_type, ordinal_type > btree_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::Experimental::UniqueToken< execution_space > token;
    std::vector< btree_type > maps( token.size() );
    auto maps_ptr = &maps;
    size_type nnz = 0;

    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::count_row_nnz",
        policy_type( 0, a_nrows ), [=]( const uint64_t row, size_type& lnnz ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          int id = token.acquire();
          auto& acc = maps_ptr->operator[]( id );
          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_idx = b_rowmap( a_entries( a_idx ) );
            auto b_end = b_rowmap( a_entries( a_idx + 1 ) + 1 );
            for ( ; b_idx != b_end; b_idx += 2 ) {
              auto old = acc[ b_entries( b_idx ) ];
              acc[ b_entries( b_idx ) ]
                  = DIVERG_MACRO_MAX( old, b_entries( b_idx + 1 ) );
            }
          }

          size_type r_nnz = 0;
          ordinal_type count = 0;
          if ( !acc.empty() ) {
            auto it = acc.begin();
            ordinal_type lo = it->first;
            ordinal_type hi = it->second;
            ++it;
            for ( ; it != acc.end(); ++it ) {
              if ( it->first <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, it->first );
                hi = DIVERG_MACRO_MAX( hi, it->second );
                continue;
              }
              r_nnz += hi - lo + 1;
              lo = it->first;
              hi = it->second;
              count += 2;
            }
            r_nnz += hi - lo + 1;
            count += 2;
            acc.clear();
          }

          token.release(id);

          c_rowmap( row + 1 ) = count;
          if ( row == 0 ) c_rowmap( 0 ) = 0;
          lnnz += r_nnz;
        }, nnz );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spgemm_symbolic::computing_row_map_c",
        policy_type( 0, a_nrows ),
        [=]( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    return nnz;
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *         (ThreadRangePolicyPartition-BTreeAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TExecGrid >
  inline void
  _range_spgemm_numeric( THandle&,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC c_rowmap,
                         TEntriesDeviceViewC& c_entries,
                         TExecGrid, ThreadRangePolicyPartition, BTreeAccumulator )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::RangePolicy< execution_space > policy_type;
    typedef std::map< ordinal_type, ordinal_type > btree_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;

    Kokkos::Experimental::UniqueToken< execution_space > token;
    std::vector< btree_type > maps( token.size() );
    auto maps_ptr = &maps;

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spgemm_numeric::compute_numeric",
        policy_type( 0, a_nrows ), [=]( const uint64_t row ) {
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          int id = token.acquire();
          auto& acc = maps_ptr->operator[]( id );
          for ( ; a_idx < a_end; a_idx += 2 ) {
            auto b_idx = b_rowmap( a_entries( a_idx ) );
            auto b_end = b_rowmap( a_entries( a_idx + 1 ) + 1 );
            for ( ; b_idx != b_end; b_idx += 2 ) {
              auto old = acc[ b_entries( b_idx ) ];
              acc[ b_entries( b_idx ) ]
                  = DIVERG_MACRO_MAX( old, b_entries( b_idx + 1 ) );
            }
          }

          if ( !acc.empty() ) {
            size_type c_idx = c_rowmap( row );
            auto it = acc.begin();
            ordinal_type lo = it->first;
            ordinal_type hi = it->second;
            ++it;
            for ( ; it != acc.end(); ++it ) {
              if ( it->first <= hi + 1 ) {  // merge
                lo = DIVERG_MACRO_MIN( lo, it->first );
                hi = DIVERG_MACRO_MAX( hi, it->second );
                continue;
              }
              c_entries( c_idx++ ) = lo;
              c_entries( c_idx++ ) = hi;
              lo = it->first;
              hi = it->second;
            }
            c_entries( c_idx++ ) = lo;
            c_entries( c_idx++ ) = hi;
            acc.clear();
          }
          token.release( id );
        } );
  }

  /**
   *  @brief Compute off-diagonal band of matrix c as the product of a and b
   *  (ThreadSequentialPartition specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: The band extremities are aligned to the bitset width.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TExecGrid, unsigned int TL1Size >
  inline typename TEntriesDeviceViewA::non_const_value_type
  rspgemm_compute_band( THandle& handle,
                        TRowMapDeviceViewA a_rowmap,
                        TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap,
                        TEntriesDeviceViewB b_entries,
                        TExecGrid grid,
                        ThreadSequentialPartition part,
                        HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewA  a_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename a_row_map_type::non_const_value_type size_type;
    typedef typename a_row_map_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto work_size = grid.team_work_size( rdensity );
    auto nof_teams = a_nrows / work_size + 1;
    auto policy = policy_type( nof_teams, team_size, vector_size );

    ordinal_type bandwidth;
    handle.init_c_min_col_index( a_nrows );
    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::compute_band", policy,
        KOKKOS_LAMBDA( const member_type& tm, ordinal_type& t_bandwidth ) {
          auto rank = tm.league_rank();
          size_type a_row = rank * work_size;
          size_type a_last_row = a_row + work_size;
          a_last_row = DIVERG_MACRO_MIN( a_last_row, a_nrows );
          ordinal_type lt_bandwidth;
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange( tm, a_row, a_last_row ),
              [ = ]( const uint64_t row, ordinal_type& r_bandwidth ) {
                auto a_idx = a_rowmap( row );
                auto a_end = a_rowmap( row + 1 );

                ordinal_type lrc_min = std::numeric_limits< ordinal_type >::max();
                ordinal_type lrc_max = std::numeric_limits< ordinal_type >::min();
                for ( ; a_idx != a_end; a_idx += 2 ) {
                  auto b_row = a_entries( a_idx );
                  auto b_last_row = a_entries( a_idx + 1 );
                  ordinal_type prc_min;
                  ordinal_type prc_max;
                  Kokkos::parallel_reduce(
                      Kokkos::ThreadVectorRange( tm, b_row, b_last_row + 1 ),
                      [=]( const uint64_t j, ordinal_type& brc_min,
                             ordinal_type& brc_max ) {
                        auto b_idx = b_rowmap( j );
                        auto b_end = b_rowmap( j + 1 );

                        if ( b_idx == b_end ) {
                          brc_min = std::numeric_limits< ordinal_type >::max();
                          brc_max = std::numeric_limits< ordinal_type >::min();
                          return;
                        }

                        auto b_min = b_entries( b_idx );
                        if ( b_min < brc_min ) {
                          brc_min = hbv_type::aligned_index( b_min );
                        }
                        auto b_max = b_entries( b_end - 1 ) + 1;
                        if ( brc_max < b_max ) {
                          brc_max = hbv_type::aligned_index_ceil( b_max );
                        }
                      },
                      Kokkos::Min< ordinal_type >( prc_min ),
                      Kokkos::Max< ordinal_type >( prc_max ) );
                  if ( prc_min < lrc_min ) lrc_min = prc_min;
                  if ( prc_max > lrc_max ) lrc_max = prc_max;
                }
                Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                  handle.c_min_col_index( row )
                      = ( lrc_min < lrc_max ) ? lrc_min : 0;
                } );
                ordinal_type lr_bandwidth = ( lrc_min < lrc_max )
                                                ? lrc_max - lrc_min
                                                : hbv_type::L1_SIZE;
                if ( lr_bandwidth > r_bandwidth ) r_bandwidth = lr_bandwidth;
              },
              Kokkos::Max< ordinal_type >( lt_bandwidth ) );
          if ( lt_bandwidth > t_bandwidth ) t_bandwidth = lt_bandwidth;
        },
        Kokkos::Max< ordinal_type >( bandwidth ) );

    assert( bandwidth >= hbv_type::L1_SIZE );
    return bandwidth;
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (ThreadSequentialPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of size
   *  `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TExecGrid,
            unsigned int TL1Size >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  _range_spgemm_symbolic( THandle& handle,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          TExecGrid grid,
                          ThreadSequentialPartition part,
                          HBitVectorAccumulator< TL1Size > acc_tag )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto work_size = grid.team_work_size( rdensity );
    auto nof_teams = a_nrows / work_size + 1;
    auto policy = policy_type( nof_teams, team_size, vector_size );

    ordinal_type bandwidth
        = rspgemm_compute_band( handle, a_rowmap, a_entries, b_rowmap,
                                b_entries, grid, part, acc_tag );
    handle.c_bandwidth = bandwidth;
    hbv_type::set_scratch_size( policy, bandwidth, part );

    size_type nnz = 0;
    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::count_row_nnz", policy,
        KOKKOS_LAMBDA( const member_type& tm, size_type& lnnz ) {
          auto rank = tm.league_rank();
          size_type a_row = rank * work_size;
          size_type a_last_row = a_row + work_size;
          a_last_row = DIVERG_MACRO_MIN( a_last_row, a_nrows );
          hbv_type hbv( tm, bandwidth, part );

          size_type t_nnz = 0;  // team nnz
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange( tm, a_row, a_last_row ),
              [ & ]( const uint64_t row, size_type& lt_nnz ) {
                auto a_idx = a_rowmap( row );
                auto a_end = a_rowmap( row + 1 );
                hbv.set_l1_at( handle.c_min_col_index( row ) );
                // min entry (bitset aligned) in the current `row` in C
                ordinal_type c_min = hbv.l1_begin_idx();
                // max entry + 1 (bitset aligned) in the current `row` in C
                ordinal_type c_max = c_min + hbv.l1_size();
                size_type r_nnz = 0;

                // Setting all L1 bitsets in `h_bv` to zero
                hbv.clear_l1( tm, part );

                for ( ; a_idx != a_end; a_idx += 2 ) {
                  auto b_row = a_entries( a_idx );
                  auto b_last_row = a_entries( a_idx + 1 );
                  for ( ; b_row <= b_last_row; ++b_row ) {
                    auto b_idx = b_rowmap( b_row );
                    auto b_end = b_rowmap( b_row + 1 );

                    if ( b_idx == b_end ) continue;

                    // Incrementally zero-initialise bitsets in L2
                    auto b_min = b_entries( b_idx );
                    // Considering the initial value of c_min, `b_min < c_min`
                    // implies that the extended range is definitely in L2.
                    if ( b_min < c_min ) {
                      b_min = hbv_type::aligned_index( b_min );
                      hbv.clear_l2_by_idx( tm, b_min, c_min, part );
                      c_min = b_min;  // update c_min
                    }
                    auto b_max = b_entries( b_end - 1 ) + 1;
                    // Considering the initial value of c_max, `c_max < b_max`
                    // implies that the extended range is definitely in L2.
                    if ( c_max < b_max ) {
                      b_max = hbv_type::aligned_index_ceil( b_max );
                      hbv.clear_l2_by_idx( tm, c_max, b_max, part );
                      c_max = b_max;  // update c_max
                    }

                    for ( ; b_idx != b_end; b_idx += 2 ) {
                      auto s = b_entries( b_idx );
                      auto f = b_entries( b_idx + 1 );
                      hbv.set( tm, s, f, part );
                    }
                  }
                }

                auto c_lbs = hbv_type::bindex( c_min );
                auto c_rbs = hbv_type::bindex( c_max );
                ordinal_type count = 0;
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange( tm, c_lbs, c_rbs ),
                    [=]( const uint64_t j, ordinal_type& local_count,
                           size_type& lr_nnz ) {
                      auto c
                          = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                      auto x = hbv( j );
                      local_count += 2 * hbv_type::cnt01( x, c );
                      lr_nnz += hbv_type::cnt( x );
                    },
                    count, r_nnz );

                Kokkos::single( Kokkos::PerThread( tm ),
                                [=, lt_nnz_ptr = &lt_nnz]() {
                                  c_rowmap( row + 1 ) = count;
                                  if ( row == 0 ) c_rowmap( 0 ) = 0;
                                  *lt_nnz_ptr += r_nnz;
                                } );
              },
              t_nnz );
          Kokkos::single( Kokkos::PerTeam( tm ),
                          [=, lnnz_ptr = &lnnz]() { *lnnz_ptr += t_nnz; } );
        },
        nnz );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spgemm_symbolic::computing_row_map_c", a_nrows,
        KOKKOS_LAMBDA( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    return nnz;
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *  (ThreadSequentialPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TExecGrid, unsigned int TL1Size >
  inline void
  _range_spgemm_numeric( THandle& handle,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC c_rowmap,
                         TEntriesDeviceViewC& c_entries,
                         TExecGrid grid,
                         ThreadSequentialPartition part,
                         HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto bandwidth = handle.c_bandwidth;
    DIVERG_ASSERT( bandwidth != 0 );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto work_size = grid.team_work_size( rdensity );
    auto nof_teams = a_nrows / work_size + 1;
    auto policy = policy_type( nof_teams, team_size, vector_size );
    hbv_type::set_scratch_size( policy, bandwidth, part );

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spgemm_numeric::accumulate_hbv", policy,
        KOKKOS_LAMBDA( const member_type& tm ) {
          auto rank = tm.league_rank();
          size_type a_row = rank * work_size;
          size_type a_last_row = a_row + work_size;
          a_last_row = DIVERG_MACRO_MIN( a_last_row, a_nrows );
          hbv_type hbv( tm, bandwidth, part );

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, a_row, a_last_row ),
              [&]( const uint64_t row ) {
                auto a_idx = a_rowmap( row );
                auto a_end = a_rowmap( row + 1 );
                hbv.set_l1_at( handle.c_min_col_index( row ) );
                // min entry (bitset aligned) in the current `row` in C
                ordinal_type c_min = hbv.l1_begin_idx();
                // max entry + 1 (bitset aligned) in the current `row` in C
                ordinal_type c_max = c_min + hbv.l1_size();

                // Setting all L1 bitsets in `h_bv` to zero
                hbv.clear_l1( tm, part );

                for ( ; a_idx != a_end; a_idx += 2 ) {
                  auto b_row = a_entries( a_idx );
                  auto b_last_row = a_entries( a_idx + 1 );
                  for ( ; b_row <= b_last_row; ++b_row ) {
                    auto b_idx = b_rowmap( b_row );
                    auto b_end = b_rowmap( b_row + 1 );

                    if ( b_idx == b_end ) continue;

                    // Incrementally zero-initialise bitsets in L2
                    auto b_min = b_entries( b_idx );
                    // Considering the initial value of c_min, `b_min < c_min`
                    // implies that the extended range is definitely in L2.
                    if ( b_min < c_min ) {
                      b_min = hbv_type::aligned_index( b_min );
                      hbv.clear_l2_by_idx( tm, b_min, c_min, part );
                      c_min = b_min;  // update c_min
                    }
                    auto b_max = b_entries( b_end - 1 ) + 1;
                    // Considering the initial value of c_max, `c_max < b_max`
                    // implies that the extended range is definitely in L2.
                    if ( c_max < b_max ) {
                      b_max = hbv_type::aligned_index_ceil( b_max );
                      hbv.clear_l2_by_idx( tm, c_max, b_max, part );
                      c_max = b_max;  // update c_max
                    }

                    for ( ; b_idx != b_end; b_idx += 2 ) {
                      auto s = b_entries( b_idx );
                      auto f = b_entries( b_idx + 1 );
                      hbv.set( tm, s, f, part );
                    }
                  }
                }

                auto c_lbs = hbv_type::bindex( c_min );
                auto c_rbs = hbv_type::bindex( c_max );
                Kokkos::parallel_for(
                    Kokkos::ThreadVectorRange( tm, c_lbs, c_rbs ),
                    [=]( const uint64_t j ) {
                      auto c
                          = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                      auto x = hbv( j );
                      auto bounds
                          = hbv_type::map01( x, c ) | hbv_type::map10( x, c );

                      if ( bounds != 0 ) {
                        size_type c_idx = 0;
                        auto c = 0;
                        for ( uint64_t k = c_lbs; k < j; ++k ) {
                          auto x = hbv( k );
                          c_idx += hbv_type::cnt01( x, c )
                                   + hbv_type::cnt10( x, c );
                          c = hbv_type::msb( hbv( k ) );
                        }
                        c_idx += c_rowmap( row );

                        for ( uint64_t k = 0; k < hbv_type::cnt( bounds ); ++k ) {
                          auto lidx = c_idx + k;
                          c_entries( lidx ) = hbv_type::start_index( j )
                                              + hbv_type::sel( bounds, k + 1 )
                                              - lidx % 2;
                        }
                      }
                    } );

                Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                  if ( hbv_type::msb( hbv( c_rbs - 1 ) ) ) {
                    c_entries( c_rowmap( row + 1 ) - 1 )
                        = hbv_type::start_index( c_rbs ) - 1;
                  }
                } );
              } );
        } );
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (TeamSequentialPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of size
   *  `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TExecGrid,
            unsigned int TL1Size >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  _range_spgemm_symbolic( THandle& handle,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          TExecGrid grid,
                          TeamSequentialPartition part,
                          HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_rnnz = b_entries.extent( 0 ) / 2;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    auto b_range_density = grid.row_density( b_row_density, b_rnnz );
    size_type bitset_count = grid.row_density( b_range_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    size_type nnz = 0;
    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::count_row_nnz", policy,
        KOKKOS_LAMBDA( const member_type& tm, size_type& lnnz ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min = hbv.l1_begin_idx();
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max = c_min + hbv.l1_size();
          size_type r_nnz = 0;

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          tm.team_barrier();

          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );
            for ( ; b_row <= b_last_row; ++b_row ) {
              auto b_idx = b_rowmap( b_row );
              auto b_end = b_rowmap( b_row + 1 );

              if ( b_idx == b_end ) continue;

              // Incrementally zero-initialise bitsets in L2
              auto b_min = b_entries( b_idx );
              // Considering the initial value of c_min, `b_min < c_min`
              // implies that the extended range is definitely in L2.
              if ( b_min < c_min ) {
                b_min = hbv_type::aligned_index( b_min );
                hbv.clear_l2_by_idx( tm, b_min, c_min, part );
                c_min = b_min;  // update c_min
              }
              auto b_max = b_entries( b_end - 1 ) + 1;
              // Considering the initial value of c_max, `c_max < b_max`
              // implies that the extended range is definitely in L2.
              if ( c_max < b_max ) {
                b_max = hbv_type::aligned_index_ceil( b_max );
                hbv.clear_l2_by_idx( tm, c_max, b_max, part );
                c_max = b_max;  // update c_max
              }

              tm.team_barrier();

              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange( tm, b_idx / 2, b_end / 2 ),
                  [&]( const uint64_t jj ) {
                    auto j = jj * 2;
                    auto s = b_entries( j );
                    auto f = b_entries( j + 1 );
                    hbv.set( tm, s, f, part );
                  } );

              tm.team_barrier();  // NOTE: necessary to avoid data race in the
                                  // next iteration
            }
          }

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          ordinal_type count = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j, ordinal_type& local_count,
                     size_type& lr_nnz ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                local_count += 2 * hbv_type::cnt01( x, c );
                lr_nnz += hbv_type::cnt( x );
              },
              count, r_nnz );

          Kokkos::single( Kokkos::PerTeam( tm ), [=, lnnz_ptr = &lnnz]() {
            c_rowmap( row + 1 ) = count;
            if ( row == 0 ) c_rowmap( 0 ) = 0;
            *lnnz_ptr += r_nnz;
          } );
        }, nnz );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spgemm_symbolic::computing_row_map_c", a_nrows,
        KOKKOS_LAMBDA( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    return nnz;
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *  (TeamSequentialPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TExecGrid, unsigned int TL1Size >
  inline void
  _range_spgemm_numeric( THandle& handle,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC c_rowmap,
                         TEntriesDeviceViewC& c_entries,
                         TExecGrid grid,
                         TeamSequentialPartition part,
                         HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_rnnz = b_entries.extent( 0 ) / 2;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    auto b_range_density = grid.row_density( b_row_density, b_rnnz );
    size_type bitset_count = grid.row_density( b_range_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spgemm_numeric::accumulate_hbv", policy,
        KOKKOS_LAMBDA( const member_type& tm ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min = hbv.l1_begin_idx();
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max = c_min + hbv.l1_size();

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          tm.team_barrier();

          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );
            for ( ; b_row <= b_last_row; ++b_row ) {
              auto b_idx = b_rowmap( b_row );
              auto b_end = b_rowmap( b_row + 1 );

              if ( b_idx == b_end ) continue;

              // Incrementally zero-initialise bitsets in L2
              auto b_min = b_entries( b_idx );
              // Considering the initial value of c_min, `b_min < c_min`
              // implies that the extended range is definitely in L2.
              if ( b_min < c_min ) {
                b_min = hbv_type::aligned_index( b_min );
                hbv.clear_l2_by_idx( tm, b_min, c_min, part );
                c_min = b_min;  // update c_min
              }
              auto b_max = b_entries( b_end - 1 ) + 1;
              // Considering the initial value of c_max, `c_max < b_max`
              // implies that the extended range is definitely in L2.
              if ( c_max < b_max ) {
                b_max = hbv_type::aligned_index_ceil( b_max );
                hbv.clear_l2_by_idx( tm, c_max, b_max, part );
                c_max = b_max;  // update c_max
              }

              tm.team_barrier();

              Kokkos::parallel_for(
                  Kokkos::TeamThreadRange( tm, b_idx / 2, b_end / 2 ),
                  [&]( const uint64_t jj ) {
                    auto j = jj * 2;
                    auto s = b_entries( j );
                    auto f = b_entries( j + 1 );
                    hbv.set( tm, s, f, part );
                  } );

              tm.team_barrier();  // NOTE: necessary to avoid data race in the
                                  // next iteration
            }
          }

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                auto bounds = hbv_type::map01( x, c ) | hbv_type::map10( x, c );

                if ( bounds != 0 ) {
                  size_type c_idx = 0;
                  Kokkos::parallel_reduce(
                      Kokkos::ThreadVectorRange( tm, c_lbs, j ),
                      [=]( const uint64_t k, size_type& lci ) {
                        auto c = ( k != c_lbs )
                                     ? hbv_type::msb( hbv( k - 1 ) )
                                     : 0;
                        auto x = hbv( k );
                        lci += hbv_type::cnt01( x, c )
                               + hbv_type::cnt10( x, c );
                      },
                      c_idx );
                  c_idx += c_rowmap( row );

                  Kokkos::parallel_for(
                      Kokkos::ThreadVectorRange( tm, hbv_type::cnt( bounds ) ),
                      [=]( const uint64_t k ) {
                        auto lidx = c_idx + k;
                        c_entries( lidx ) = hbv_type::start_index( j )
                                            + hbv_type::sel( bounds, k + 1 )
                                            - lidx % 2;
                      } );
                }
              } );

          Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
            if ( hbv_type::msb( hbv( c_rbs - 1 ) ) ) {
              c_entries( c_rowmap( row + 1 ) - 1 )
                  = hbv_type::start_index( c_rbs ) - 1;
            }
          } );
        } );
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (ThreadParallelPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of size
   *  `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TExecGrid,
            unsigned int TL1Size >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  _range_spgemm_symbolic( THandle& handle,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          TExecGrid grid,
                          ThreadParallelPartition part,
                          HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    size_type nnz = 0;
    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::count_row_nnz", policy,
        KOKKOS_LAMBDA( const member_type& tm, size_type& lnnz ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          ordinal_type l1b = hbv.l1_begin_idx();
          ordinal_type l1e = l1b + hbv.l1_size();
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min;
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max;

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          // Calculating min and max column indices in the 'row'
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange( tm, a_idx / 2, a_end / 2 ),
              [&]( const uint64_t ii,
                   ordinal_type& lc_min, ordinal_type& lc_max ) {
                auto i = ii * 2;
                auto b_row = a_entries( i );
                auto b_last_row = a_entries( i + 1 );
                ordinal_type br_min;  // B row range min
                ordinal_type br_max;  // B row range max
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange( tm, b_row, b_last_row + 1 ),
                    [&]( const uint64_t j,
                         ordinal_type& lbr_min, ordinal_type& lbr_max ) {
                      auto b_idx = b_rowmap( j );
                      auto b_end = b_rowmap( j + 1 );

                      if ( b_idx == b_end ) return;

                      auto b_min = b_entries( b_idx );
                      if ( b_min < lbr_min ) {
                        b_min = hbv_type::aligned_index( b_min );
                        lbr_min = b_min;  // update lbr_min
                      }

                      auto b_max = b_entries( b_end - 1 ) + 1;
                      if ( lbr_max < b_max ) {
                        b_max = hbv_type::aligned_index_ceil( b_max );
                        lbr_max = b_max;  // update br_max
                      }
                    },
                    Kokkos::Min< ordinal_type >( br_min ),
                    Kokkos::Max< ordinal_type >( br_max ) );
                if ( br_min < lc_min ) lc_min = br_min;
                if ( br_max > lc_max ) lc_max = br_max;
              },
              Kokkos::Min< ordinal_type >( c_min ),
              Kokkos::Max< ordinal_type >( c_max ) );

          if ( c_max < c_min ) {
            c_min = l1b;
            c_max = l1b;
          }
          if ( c_min < l1b )
            hbv.clear_l2_by_idx( tm, c_min, l1b, part );
          if ( l1e < c_max )
            hbv.clear_l2_by_idx( tm, l1e, c_max, part );

          tm.team_barrier();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, a_idx / 2, a_end / 2 ),
              [&]( const uint64_t ii ) {
                auto i = ii * 2;
                auto b_row = a_entries( i );
                auto b_last_row = a_entries( i + 1 );
                auto b_idx = b_rowmap( b_row );
                auto b_end = b_rowmap( b_last_row + 1 );
                for ( ; b_idx < b_end; b_idx += 2 ) {
                  auto s = b_entries( b_idx );
                  auto f = b_entries( b_idx + 1 );
                  hbv.set( tm, s, f, part );
                }
              } );

          tm.team_barrier();

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          ordinal_type count = 0;
          size_type r_nnz = 0;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j, ordinal_type& local_count,
                   size_type& lr_nnz ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                local_count += 2 * hbv_type::cnt01( x, c );
                lr_nnz += hbv_type::cnt( x );
              },
              count, r_nnz );

          Kokkos::single( Kokkos::PerTeam( tm ), [=, lnnz_ptr = &lnnz]() {
            c_rowmap( row + 1 ) = count;
            if ( row == 0 ) c_rowmap( 0 ) = 0;
            *lnnz_ptr += r_nnz;
          } );
        }, nnz );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spgemm_symbolic::computing_row_map_c", a_nrows,
        KOKKOS_LAMBDA( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    return nnz;
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *  (ThreadParallelPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TExecGrid, unsigned int TL1Size >
  inline void
  _range_spgemm_numeric( THandle& handle,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC c_rowmap,
                         TEntriesDeviceViewC& c_entries,
                         TExecGrid grid,
                         ThreadParallelPartition part,
                         HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spgemm_numeric::accumulate_hbv", policy,
        KOKKOS_LAMBDA( const member_type& tm ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          ordinal_type l1b = hbv.l1_begin_idx();
          ordinal_type l1e = l1b + hbv.l1_size();
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min;
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max;

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          // Calculating min and max column indices in the 'row'
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange( tm, a_idx / 2, a_end / 2 ),
              [&]( const uint64_t ii,
                   ordinal_type& lc_min, ordinal_type& lc_max ) {
                auto i = ii * 2;
                auto b_row = a_entries( i );
                auto b_last_row = a_entries( i + 1 );
                ordinal_type br_min;  // B row range min
                ordinal_type br_max;  // B row range max
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange( tm, b_row, b_last_row + 1 ),
                    [&]( const uint64_t j,
                         ordinal_type& lbr_min, ordinal_type& lbr_max ) {
                      auto b_idx = b_rowmap( j );
                      auto b_end = b_rowmap( j + 1 );

                      if ( b_idx == b_end ) return;

                      auto b_min = b_entries( b_idx );
                      if ( b_min < lbr_min ) {
                        b_min = hbv_type::aligned_index( b_min );
                        lbr_min = b_min;  // update lbr_min
                      }

                      auto b_max = b_entries( b_end - 1 ) + 1;
                      if ( lbr_max < b_max ) {
                        b_max = hbv_type::aligned_index_ceil( b_max );
                        lbr_max = b_max;  // update br_max
                      }
                    },
                    Kokkos::Min< ordinal_type >( br_min ),
                    Kokkos::Max< ordinal_type >( br_max ) );
                if ( br_min < lc_min ) lc_min = br_min;
                if ( br_max > lc_max ) lc_max = br_max;
              },
              Kokkos::Min< ordinal_type >( c_min ),
              Kokkos::Max< ordinal_type >( c_max ) );

          if ( c_max < c_min ) {
            c_min = l1b;
            c_max = l1b;
          }
          if ( c_min < l1b )
            hbv.clear_l2_by_idx( tm, c_min, l1b, part );
          if ( l1e < c_max )
            hbv.clear_l2_by_idx( tm, l1e, c_max, part );

          tm.team_barrier();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, a_idx / 2, a_end / 2 ),
              [&]( const uint64_t ii ) {
                auto i = ii * 2;
                auto b_row = a_entries( i );
                auto b_last_row = a_entries( i + 1 );
                auto b_idx = b_rowmap( b_row );
                auto b_end = b_rowmap( b_last_row + 1 );
                for ( ; b_idx < b_end; b_idx += 2 ) {
                  auto s = b_entries( b_idx );
                  auto f = b_entries( b_idx + 1 );
                  hbv.set( tm, s, f, part );
                }
              } );

          tm.team_barrier();

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                auto bounds = hbv_type::map01( x, c ) | hbv_type::map10( x, c );

                if ( bounds != 0 ) {
                  size_type c_idx = 0;
                  Kokkos::parallel_reduce(
                      Kokkos::ThreadVectorRange( tm, c_lbs, j ),
                      [=]( const uint64_t k, size_type& lci ) {
                        auto c = ( k != c_lbs )
                                     ? hbv_type::msb( hbv( k - 1 ) )
                                     : 0;
                        auto x = hbv( k );
                        lci += hbv_type::cnt01( x, c )
                               + hbv_type::cnt10( x, c );
                      },
                      c_idx );
                  c_idx += c_rowmap( row );

                  Kokkos::parallel_for(
                      Kokkos::ThreadVectorRange( tm, hbv_type::cnt( bounds ) ),
                      [=]( const uint64_t k ) {
                        auto lidx = c_idx + k;
                        c_entries( lidx ) = hbv_type::start_index( j )
                                            + hbv_type::sel( bounds, k + 1 )
                                            - lidx % 2;
                      } );
                }
              } );

          Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
            if ( c_lbs < c_rbs && hbv_type::msb( hbv( c_rbs - 1 ) ) ) {
              c_entries( c_rowmap( row + 1 ) - 1 )
                  = hbv_type::start_index( c_rbs ) - 1;
            }
          } );
        } );
  }

  /**
   *  @brief Symbolic phase of computing matrix c as the product of a and b
   *  (ThreadRangeParallelPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` is allocated on device and is of size
   *  `a.numRows()+1`.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TExecGrid,
            unsigned int TL1Size >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  _range_spgemm_symbolic( THandle& handle,
                          TRowMapDeviceViewA a_rowmap,
                          TEntriesDeviceViewA a_entries,
                          TRowMapDeviceViewB b_rowmap,
                          TEntriesDeviceViewB b_entries,
                          TRowMapDeviceViewC& c_rowmap,
                          TExecGrid grid,
                          ThreadRangeParallelPartition part,
                          HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename a_entries_type::non_const_value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_row_map_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    size_type nnz = 0;
    Kokkos::parallel_reduce(
        "diverg::crs_matrix::range_spgemm_symbolic::count_row_nnz", policy,
        KOKKOS_LAMBDA( const member_type& tm, size_type& lnnz ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          ordinal_type l1b = hbv.l1_begin_idx();
          ordinal_type l1e = l1b + hbv.l1_size();
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min = l1b;
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max = l1e;

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          // Calculating min and max column indices in the 'row'
          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );

            ordinal_type br_min;  // B row range min
            ordinal_type br_max;  // B row range max

            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange( tm, b_row, b_last_row + 1 ),
                [&]( const uint64_t i,
                     ordinal_type& lbr_min, ordinal_type& lbr_max ) {
                  auto b_idx = b_rowmap( i );
                  auto b_end = b_rowmap( i + 1 );

                  if ( b_idx == b_end ) return;

                  auto b_min = b_entries( b_idx );
                  if ( b_min < lbr_min ) {
                    b_min = hbv_type::aligned_index( b_min );
                    lbr_min = b_min;  // update lbr_min
                  }
                  auto b_max = b_entries( b_end - 1 ) + 1;
                  if ( lbr_max < b_max ) {
                    b_max = hbv_type::aligned_index_ceil( b_max );
                    lbr_max = b_max;  // update lbr_max
                  }
                },
                Kokkos::Min< ordinal_type >( br_min ),
                Kokkos::Max< ordinal_type >( br_max ) );

            if ( br_min < c_min ) c_min = br_min;
            if ( br_max > c_max ) c_max = br_max;
          }

          if ( c_min < l1b )
            hbv.clear_l2_by_idx( tm, c_min, l1b, part );
          if ( l1e < c_max )
            hbv.clear_l2_by_idx( tm, l1e, c_max, part );

          tm.team_barrier();

          for ( a_idx = a_rowmap( row ); a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );
            auto b_idx = b_rowmap( b_row );
            auto b_end = b_rowmap( b_last_row + 1 );
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, b_idx / 2, b_end / 2 ),
                [&]( const uint64_t ii ) {
                  auto i = ii * 2;
                  auto s = b_entries( i );
                  auto f = b_entries( i + 1 );
                  hbv.set( tm, s, f, part );
                } );
          }

          tm.team_barrier();

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          ordinal_type count = 0;
          size_type r_nnz;
          Kokkos::parallel_reduce(
              Kokkos::TeamVectorRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j, ordinal_type& local_count,
                   size_type& lr_nnz ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                local_count += 2 * hbv_type::cnt01( x, c );
                lr_nnz += hbv_type::cnt( x );
              },
              count, r_nnz );

          Kokkos::single( Kokkos::PerTeam( tm ), [=, lnnz_ptr = &lnnz]() {
            c_rowmap( row + 1 ) = count;
            if ( row == 0 ) c_rowmap( 0 ) = 0;
            *lnnz_ptr += r_nnz;
          } );
        }, nnz );

    Kokkos::parallel_scan(
        "diverg::crs_matrix::range_spgemm_symbolic::computing_row_map_c", a_nrows,
        KOKKOS_LAMBDA( const int i, size_type& update, const bool final ) {
          // Load old value in case we update it before accumulating
          const size_type val_ip1 = c_rowmap( i + 1 );
          update += val_ip1;
          if ( final )
            c_rowmap( i + 1 ) = update;  // only update array on final pass
        } );

    return nnz;
  }

  /**
   *  @brief Numeric phase of computing matrix c as the product of a and b
   *  (ThreadRangeParallelPartition-HBitVectorAccumulator specialisation).
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   *  NOTE: This function assumes `c_rowmap` and `c_entries` are allocated on
   *        device with sufficient space.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TExecGrid, unsigned int TL1Size >
  inline void
  _range_spgemm_numeric( THandle& handle,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC c_rowmap,
                         TEntriesDeviceViewC& c_entries,
                         TExecGrid grid,
                         ThreadRangeParallelPartition part,
                         HBitVectorAccumulator< TL1Size > )
  {
    typedef TEntriesDeviceViewA a_entries_type;
    typedef TEntriesDeviceViewB b_entries_type;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;
    typedef typename c_entries_type::execution_space execution_space;
    typedef Kokkos::TeamPolicy< execution_space > policy_type;
    typedef typename policy_type::member_type member_type;
    typedef diverg::HBitVector< TL1Size, execution_space > hbv_type;

    // TODO: Extend static asserts to all views
    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename b_entries_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    static_assert(
        std::is_same< typename a_entries_type::memory_space,
                      typename c_row_map_type::memory_space >::value,
        "both entries and row map views should be in the same memory space" );

    auto a_nrows = a_rowmap.extent( 0 ) - 1;
    auto b_nrows = b_rowmap.extent( 0 ) - 1;
    auto b_ncols = handle.b_ncols;
    auto b_row_density = grid.row_density( handle.b_nnz, b_nrows );
    size_type bitset_count = grid.row_density( b_row_density, hbv_type::BITSET_WIDTH );
    auto rdensity = DIVERG_MACRO_MAX( bitset_count, hbv_type::l1_num_bitsets() );

    auto vector_size = grid.vector_size( rdensity );
    auto team_size = grid.team_size( rdensity );
    auto policy = policy_type( a_nrows, team_size, vector_size );
    hbv_type::set_scratch_size( policy, b_ncols, part );

    Kokkos::parallel_for(
        "diverg::crs_matrix::range_spgemm_numeric::accumulate_hbv", policy,
        KOKKOS_LAMBDA( const member_type& tm ) {
          auto row = tm.league_rank();
          auto a_idx = a_rowmap( row );
          auto a_end = a_rowmap( row + 1 );
          hbv_type hbv( tm, b_ncols, row, part );
          ordinal_type l1b = hbv.l1_begin_idx();
          ordinal_type l1e = l1b + hbv.l1_size();
          // min entry (bitset aligned) in the current `row` in C
          ordinal_type c_min = l1b;
          // max entry + 1 (bitset aligned) in the current `row` in C
          ordinal_type c_max = l1e;

          // Setting all L1 bitsets in `h_bv` to zero
          hbv.clear_l1( tm, part );

          // Calculating min and max column indices in the 'row'
          for ( ; a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );

            ordinal_type br_min;  // B row range min
            ordinal_type br_max;  // B row range max

            Kokkos::parallel_reduce(
                Kokkos::TeamVectorRange( tm, b_row, b_last_row + 1 ),
                [&]( const uint64_t i,
                     ordinal_type& lbr_min, ordinal_type& lbr_max ) {
                  auto b_idx = b_rowmap( i );
                  auto b_end = b_rowmap( i + 1 );

                  if ( b_idx == b_end ) return;

                  auto b_min = b_entries( b_idx );
                  if ( b_min < lbr_min ) {
                    b_min = hbv_type::aligned_index( b_min );
                    lbr_min = b_min;  // update lbr_min
                  }
                  auto b_max = b_entries( b_end - 1 ) + 1;
                  if ( lbr_max < b_max ) {
                    b_max = hbv_type::aligned_index_ceil( b_max );
                    lbr_max = b_max;  // update lbr_max
                  }
                },
                Kokkos::Min< ordinal_type >( br_min ),
                Kokkos::Max< ordinal_type >( br_max ) );

            if ( br_min < c_min ) c_min = br_min;
            if ( br_max > c_max ) c_max = br_max;
          }

          if ( c_min < l1b )
            hbv.clear_l2_by_idx( tm, c_min, l1b, part );
          if ( l1e < c_max )
            hbv.clear_l2_by_idx( tm, l1e, c_max, part );

          tm.team_barrier();

          for ( a_idx = a_rowmap( row ); a_idx != a_end; a_idx += 2 ) {
            auto b_row = a_entries( a_idx );
            auto b_last_row = a_entries( a_idx + 1 );
            auto b_idx = b_rowmap( b_row );
            auto b_end = b_rowmap( b_last_row + 1 );
            Kokkos::parallel_for(
                Kokkos::TeamThreadRange( tm, b_idx / 2, b_end / 2 ),
                [&]( const uint64_t ii ) {
                  auto i = ii * 2;
                  auto s = b_entries( i );
                  auto f = b_entries( i + 1 );
                  hbv.set( tm, s, f, part );
                } );
          }

          tm.team_barrier();

          auto c_lbs = hbv_type::bindex( c_min );
          auto c_rbs = hbv_type::bindex( c_max );
          Kokkos::parallel_for(
              Kokkos::TeamThreadRange( tm, c_lbs, c_rbs ),
              [=]( const uint64_t j ) {
                auto c = ( j != c_lbs ) ? hbv_type::msb( hbv( j - 1 ) ) : 0;
                auto x = hbv( j );
                auto bounds = hbv_type::map01( x, c ) | hbv_type::map10( x, c );

                if ( bounds != 0 ) {
                  size_type c_idx = 0;
                  Kokkos::parallel_reduce(
                      Kokkos::ThreadVectorRange( tm, c_lbs, j ),
                      [=]( const uint64_t k, size_type& lci ) {
                        auto c = ( k != c_lbs )
                                     ? hbv_type::msb( hbv( k - 1 ) )
                                     : 0;
                        auto x = hbv( k );
                        lci += hbv_type::cnt01( x, c )
                               + hbv_type::cnt10( x, c );
                      },
                      c_idx );
                  c_idx += c_rowmap( row );

                  Kokkos::parallel_for(
                      Kokkos::ThreadVectorRange( tm, hbv_type::cnt( bounds ) ),
                      [=]( const uint64_t k ) {
                        auto lidx = c_idx + k;
                        c_entries( lidx ) = hbv_type::start_index( j )
                                            + hbv_type::sel( bounds, k + 1 )
                                            - lidx % 2;
                      } );
                }
              } );

          Kokkos::single( Kokkos::PerTeam( tm ), [=]() {
            if ( c_lbs < c_rbs && hbv_type::msb( hbv( c_rbs - 1 ) ) ) {
              c_entries( c_rowmap( row + 1 ) - 1 )
                  = hbv_type::start_index( c_rbs ) - 1;
            }
          } );
        } );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  range_spgemm_symbolic( THandle& handle,
                         TRowMapDeviceViewA a_rowmap,
                         TEntriesDeviceViewA a_entries,
                         TRowMapDeviceViewB b_rowmap,
                         TEntriesDeviceViewB b_entries,
                         TRowMapDeviceViewC& c_rowmap,
                         TSparseConfig config={} )
  {
    typedef typename TEntriesDeviceViewA::non_const_value_type ordinal_type;
    //typedef typename TSparseConfig::partition_type partition_type;
    //typedef typename TSparseConfig::accumulator_type accumulator_type;
    //typedef typename TSparseConfig::grid_type grid_type;

    assert( handle.a_ncols == static_cast< ordinal_type >( b_rowmap.extent( 0 ) - 1 ) );
    DIVERG_ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    DIVERG_ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    return _range_spgemm_symbolic( handle, a_rowmap, a_entries, b_rowmap,
                                   b_entries, c_rowmap, config.grid,
                                   config.part, config.accm );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline void
  range_spgemm_numeric( THandle& handle,
                        TRowMapDeviceViewA a_rowmap, TEntriesDeviceViewA a_entries,
                        TRowMapDeviceViewB b_rowmap, TEntriesDeviceViewB b_entries,
                        TRowMapDeviceViewC c_rowmap, TEntriesDeviceViewC& c_entries,
                        TSparseConfig config={} )
  {
    typedef typename TEntriesDeviceViewC::value_type ordinal_type;
    //typedef typename TSparseConfig::partition_type partition_type;
    //typedef typename TSparseConfig::accumulator_type accumulator_type;
    //typedef typename TSparseConfig::grid_type grid_type;

    assert( handle.a_ncols == static_cast< ordinal_type >( b_rowmap.extent( 0 ) - 1 ) );
    DIVERG_ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    DIVERG_ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    _range_spgemm_numeric( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                           c_rowmap, c_entries, config.grid, config.part,
                           config.accm );
  }

  /**
   *  @brief Computing matrix c as the product of a and b.
   *
   *  NOTE: All matrices are assumed to be in Range CRS format.
   */
  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewB, typename TEntriesDeviceViewB,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TSparseConfig=DefaultSparseConfiguration,
            typename TTimer = Kokkos::Timer >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  range_spgemm( THandle& handle,
                TRowMapDeviceViewA a_rowmap, TEntriesDeviceViewA a_entries,
                TRowMapDeviceViewB b_rowmap, TEntriesDeviceViewB b_entries,
                TRowMapDeviceViewC& c_rowmap, TEntriesDeviceViewC& c_entries,
                TSparseConfig config={},
                TTimer* timer_ptr = nullptr )
  {
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_entries_type::value_type ordinal_type;
    typedef typename c_row_map_type::value_type size_type;

    assert( handle.a_ncols == static_cast< ordinal_type >( b_rowmap.extent( 0 ) - 1 ) );
    DIVERG_ASSERT( handle.a_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );
    DIVERG_ASSERT( handle.b_ncols <= std::numeric_limits< ordinal_type >::max() - 1 );

    ordinal_type n = a_rowmap.extent( 0 ) - 1;
    c_rowmap = c_row_map_type(
        Kokkos::ViewAllocateWithoutInitializing( "c_rowmap" ),
        a_rowmap.extent( 0 ) );

#ifdef DIVERG_STATS
    if ( timer_ptr ) timer_ptr->reset();
#endif

    auto nnz = range_spgemm_symbolic( handle, a_rowmap, a_entries, b_rowmap,
                                      b_entries, c_rowmap, config );

#ifdef DIVERG_STATS
    if ( timer_ptr ) {
      config.space.fence();
      auto d = timer_ptr->seconds();
      std::cout << "diverg::range_spgemm_symbolic time: " << d * 1000 << "ms"
                << std::endl;
    }
#endif

    size_type c_rnnz;
    Kokkos::deep_copy( c_rnnz, Kokkos::subview( c_rowmap, n ) );
    c_entries = c_entries_type( Kokkos::ViewAllocateWithoutInitializing( "C" ),
                                c_rnnz );

#ifdef DIVERG_STATS
    if ( timer_ptr ) timer_ptr->reset();
#endif

    range_spgemm_numeric( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                          c_rowmap, c_entries, config );

#ifdef DIVERG_STATS
    if ( timer_ptr ) {
      config.space.fence();
      auto d = timer_ptr->seconds();
      std::cout << "diverg::range_spgemm_numeric time: " << d * 1000 << "ms"
                << std::endl;
    }
#endif

    return nnz;
  }

  template< typename TRCRSMatrix,
            typename TSparseConfig=DefaultSparseConfiguration >
  inline TRCRSMatrix
  range_spgemm( TRCRSMatrix const& a, TRCRSMatrix const& b,
                TSparseConfig config={} )
  {
    typedef TRCRSMatrix range_crsmatrix_t;
    typedef typename range_crsmatrix_t::ordinal_type ordinal_type;
    //typedef TSparseConfig config_type;
    //typedef typename config_type::execution_space execution_space;

    assert( a.numCols() == b.numRows() );
    DIVERG_ASSERT( a.numCols() <= std::numeric_limits< ordinal_type >::max() - 1 );
    DIVERG_ASSERT( b.numCols() <= std::numeric_limits< ordinal_type >::max() - 1 );

    auto a_entries = a.entries_device_view( config.space );
    auto a_rowmap = a.rowmap_device_view( config.space );
    auto b_entries = b.entries_device_view( config.space );
    auto b_rowmap = b.rowmap_device_view( config.space );

    auto c_entries = range_crsmatrix_t::make_entries_device_view( config.space );
    auto c_rowmap = range_crsmatrix_t::make_rowmap_device_view( config.space );

    SparseRangeHandle handle( a, b, config.space );

    auto nnz = range_spgemm( handle, a_rowmap, a_entries, b_rowmap, b_entries,
                             c_rowmap, c_entries, config );

    // FIXME: since entries and rowmap arrays of range CRS is not a view, there
    // would be an extra copy here and the `c_entries` and `c_rowmap` cannot be
    // moved when the views are on the same memory space (the RCRS ctor does call
    // `deep_copy`).
    return TRCRSMatrix( b.numCols(), c_entries, c_rowmap, nnz );
  }

  template< typename TRCRSMatrix,
            typename TSparseConfig = DefaultSparseConfiguration >
  inline TRCRSMatrix
  range_power( TRCRSMatrix const& a, unsigned int k, TSparseConfig config = {} )
  {
    typedef TRCRSMatrix rcrsmatrix_t;

    auto a_entries = a.entries_device_view( config.space );
    auto a_rowmap = a.rowmap_device_view( config.space );

    auto c_entries = rcrsmatrix_t::make_entries_device_view( config.space );
    auto c_rowmap = rcrsmatrix_t::make_rowmap_device_view( config.space );

    SparseRangeHandle handle( a, a, config.space );
    auto nnz = range_power_inplace( handle, a_rowmap, a_entries, c_rowmap,
                                    c_entries, k, config );

    return TRCRSMatrix( a.numCols(), c_entries, c_rowmap, nnz );
  }

  template< typename THandle,
            typename TRowMapDeviceViewA, typename TEntriesDeviceViewA,
            typename TRowMapDeviceViewC, typename TEntriesDeviceViewC,
            typename TSparseConfig = DefaultSparseConfiguration,
            typename TTimer = Kokkos::Timer >
  inline typename TRowMapDeviceViewC::value_type /* size_type */
  range_power_inplace( THandle& h1,
                       TRowMapDeviceViewA& a_rowmap, TEntriesDeviceViewA& a_entries,
                       TRowMapDeviceViewC& c_rowmap, TEntriesDeviceViewC& c_entries,
                       unsigned int k, TSparseConfig config={},
                       TTimer* timer_ptr = nullptr )
  {
    typedef THandle handle_t;
    typedef TEntriesDeviceViewC c_entries_type;
    typedef TRowMapDeviceViewC  c_row_map_type;
    typedef typename c_row_map_type::value_type size_type;

    c_entries_type tmp_entries;
    c_row_map_type tmp_rowmap;

    auto ncols = h1.a_ncols;
    size_type a_nnz = h1.a_nnz;
    size_type c_nnz = ncols;

#ifdef DIVERG_STATS
    if ( timer_ptr ) timer_ptr->reset();
#endif

    create_range_identity_matrix( c_rowmap, c_entries, ncols /*== nrows*/ );

    while ( true ) {
      if ( k & 1 ) {
        handle_t handle( ncols, c_nnz, ncols, a_nnz );
        c_nnz = range_spgemm( handle, c_rowmap, c_entries, a_rowmap, a_entries,
                              tmp_rowmap, tmp_entries, config );
        c_entries = tmp_entries;
        c_rowmap = tmp_rowmap;
      }

      k = k >> 1;

      if ( k != 0 ) {
        handle_t handle( ncols, a_nnz, ncols, a_nnz );
        a_nnz = range_spgemm( handle, a_rowmap, a_entries, a_rowmap, a_entries,
                              tmp_rowmap, tmp_entries, config );
        a_entries = tmp_entries;
        a_rowmap = tmp_rowmap;
      }
      else break;
    }

#ifdef DIVERG_STATS
    if ( timer_ptr ) {
      config.space.fence();
      auto duration = timer_ptr->seconds();
      std::cout << "range_power time: " << duration * 1000 << "ms"
                << std::endl;
    }
#endif

    return c_nnz;
  }
}  // namespace diverg

#endif  // DIVERG_RANGE_SPARSE_HPP_
