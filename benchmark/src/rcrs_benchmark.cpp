/**
 *    @file  rcrs_benchmark.cpp
 *   @brief  Benchmark Range CRS matrix operations
 *
 *  Performance benchmarking of Range CRS matrix operations compared to general
 *  matrix representation as well as generic compressed row storage format (CRS).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu May 11, 2023  18:34
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include <gum/graph.hpp>
#include <gum/io_utils.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <KokkosSparse_spgemm.hpp>
#include <KokkosSparse_spadd.hpp>

#include "rcrs_benchmark.hpp"

using namespace diverg;


template< typename TXCRSMatrix1, typename TXCRSMatrix2 >
TXCRSMatrix1
copy_xcrs( TXCRSMatrix2 mat, Kokkos::Timer* timer_ptr = nullptr )
{
  using dst_matrix_type  = TXCRSMatrix1;
  using src_matrix_type  = TXCRSMatrix2;

  using ordinal_t = typename src_matrix_type::ordinal_type;
  using size_type = typename src_matrix_type::size_type;
  using scalar_t = typename src_matrix_type::value_type;

  static_assert(
      std::is_same< ordinal_t, typename dst_matrix_type::ordinal_type >::value,
      "ordinal types of two matrices should be identical" );
  static_assert(
      std::is_same< size_type, typename dst_matrix_type::size_type >::value,
      "size types of two matrices should be identical" );
  static_assert(
      std::is_same< scalar_t, typename dst_matrix_type::value_type >::value,
      "value types of two matrices should be identical" );

  using dst_exec_space = typename dst_matrix_type::execution_space;
  using dst_graph_type   = typename dst_matrix_type::staticcrsgraph_type;

  if ( timer_ptr ) timer_ptr->reset();

  auto row_map =
      Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                           mat.graph.row_map );
  auto entries =
      Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                           mat.graph.entries );
  auto values
      = Kokkos::create_mirror_view_and_copy( dst_exec_space{},
                                             mat.values );

  // using dst_row_map_type = typename dst_graph_type::row_map_type;
  // using dst_entries_type = typename dst_graph_type::entries_type;
  // using dst_values_type = typename dst_matrix_type::values_type;

  // typename dst_row_map_type::non_const_type row_map( "rowmap", mat.numRows() + 1 );
  // typename dst_entries_type::non_const_type entries( "entries", mat.nnz() );
  // typename dst_values_type::non_const_type values( "values", mat.nnz() );

  // Kokkos::deep_copy( row_map, mat.graph.row_map );
  // Kokkos::deep_copy( entries, mat.graph.entries );
  // Kokkos::deep_copy( values, mat.values );

  dst_graph_type crs_graph( entries, row_map );

  if ( timer_ptr ) {
    dst_exec_space{}.fence();
    auto duration = timer_ptr->seconds();
    std::cout << "diverg::rcrs_benchmark::copy_xcrs time: " << duration * 1000
              << "ms" << std::endl;
  }

  return dst_matrix_type( "moved", mat.numRows(), values, crs_graph );
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
kokkos_kernels_spgemm( TXCRSMatrix const& a, TXCRSMatrix const& b,
                       Kokkos::Timer* timer_ptr = nullptr )
{
  typedef typename TXCRSMatrix::ordinal_type ordinal_t;
  typedef typename TXCRSMatrix::size_type size_type;
  typedef typename TXCRSMatrix::value_type scalar_t;
  typedef typename TXCRSMatrix::execution_space execution_space;
  typedef typename TXCRSMatrix::memory_space memory_space;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
    size_type, ordinal_t, scalar_t,
    execution_space, memory_space, memory_space > kernel_handle_t;

  kernel_handle_t handle;
  handle.set_team_work_size( 16 );
  handle.set_dynamic_scheduling( true );

  KokkosSparse::SPGEMMAlgorithm spgemm_algorithm =
    // spgemm algorithm, limited by configuration at compile-time and set via the handle
    // Other options: {SPGEMM_KK_SPEED, SPGEMM_KK_MEMSPEED, SPGEMM_MKL}
      KokkosSparse::SPGEMM_KK_MEMORY;
  handle.create_spgemm_handle( spgemm_algorithm );

  TXCRSMatrix c;

  {
    if ( timer_ptr ) timer_ptr->reset();

    KokkosSparse::spgemm_symbolic( handle, a, false, b, false, c );

    if ( timer_ptr ) {
      execution_space{}.fence();
      auto duration = timer_ptr->seconds();
      std::cout << "Kokkos::SpGEMM_symbolic time: " << duration * 1000 << "ms"
                << std::endl;
    }
  }

  {
    if ( timer_ptr ) timer_ptr->reset();

    KokkosSparse::spgemm_numeric( handle, a, false, b, false, c );

    if ( timer_ptr ) {
      execution_space{}.fence();
      auto duration = timer_ptr->seconds();
      std::cout << "Kokkos::SpGEMM_numeric time: " << duration * 1000 << "ms"
                << std::endl;
    }
  }

  handle.destroy_spgemm_handle();

//  Kokkos::parallel_for(
//      "diverg::test_range_sparse::set_values",
//      Kokkos::RangePolicy< execution_space >( 0, c.nnz() ),
//      KOKKOS_LAMBDA ( const uint64_t i ) {
//        c.values( i ) = 1;
//      } );

  return c;
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
kokkos_kernels_spadd( TXCRSMatrix const& a, TXCRSMatrix const& b,
                      Kokkos::Timer* timer_ptr = nullptr )
{
  typedef typename TXCRSMatrix::ordinal_type ordinal_t;
  typedef typename TXCRSMatrix::size_type size_type;
  typedef typename TXCRSMatrix::value_type scalar_t;
  typedef typename TXCRSMatrix::execution_space execution_space;
  typedef typename TXCRSMatrix::memory_space memory_space;
  typedef typename KokkosKernels::Experimental::KokkosKernelsHandle<
    size_type, ordinal_t, scalar_t,
    execution_space, memory_space, memory_space > kernel_handle_t;

  kernel_handle_t handle;
  handle.create_spadd_handle( true /* sorted rows */ );

  TXCRSMatrix c;

  {
    if ( timer_ptr ) timer_ptr->reset();

    KokkosSparse::spadd_symbolic( &handle, a, b, c );

    if ( timer_ptr ) {
      execution_space{}.fence();
      auto duration = timer_ptr->seconds();
      std::cout << "Kokkos::SpAdd_symbolic time: " << duration * 1000 << "ms"
                << std::endl;
    }
  }

  {
    if ( timer_ptr ) timer_ptr->reset();

    KokkosSparse::spadd_numeric( &handle, 1, a, 1, b, c );

    if ( timer_ptr ) {
      execution_space{}.fence();
      auto duration = timer_ptr->seconds();
      std::cout << "Kokkos::SpAdd_numeric time: " << duration * 1000 << "ms"
                << std::endl;
    }
  }

  handle.destroy_spadd_handle();

//  Kokkos::parallel_for(
//      "diverg::test_range_sparse::set_values",
//      Kokkos::RangePolicy< execution_space >( 0, c.nnz() ),
//      KOKKOS_LAMBDA ( const uint64_t i ) {
//        c.values( i ) = 1;
//      } );

  return c;
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
kokkos_kernels_power( TXCRSMatrix const& a, unsigned int n,
                      Kokkos::Timer* timer_ptr = nullptr )
{
  using execution_space = typename TXCRSMatrix::execution_space;

  assert( a.numRows() == a.numCols() );

  if ( timer_ptr ) timer_ptr->reset();
  auto c = create_identity_matrix< TXCRSMatrix >( a.numRows() );
  auto a_copy = a;

  while ( true ) {
    if ( n & 1 ) c = kokkos_kernels_spgemm( c, a_copy );
    n = n >> 1;
    if ( n == 0 ) break;
    a_copy = kokkos_kernels_spgemm( a_copy, a_copy );
  }

  if ( timer_ptr ) {
    execution_space{}.fence();
    auto duration = timer_ptr->seconds();
    std::cout << "KokkosKernels::power time: " << duration * 1000 << "ms"
              << std::endl;
  }

  return c;
}

template< typename TXCRSMatrix >
inline TXCRSMatrix
create_dindex_pairg( TXCRSMatrix a, int dlo, int dup, int verbose = 0,
                     Kokkos::Timer* timer1_ptr = nullptr,
                     Kokkos::Timer* timer2_ptr = nullptr )
{
  using xcrsmatrix_t = TXCRSMatrix;
  using execution_space = typename xcrsmatrix_t::execution_space;

  execution_space space;

  xcrsmatrix_t c;
  {
    xcrsmatrix_t aid;
    xcrsmatrix_t ad;
    {
      xcrsmatrix_t ai;
      {
        {
          if ( timer1_ptr ) timer1_ptr->reset();
          if ( timer2_ptr ) timer2_ptr->reset();

          auto I = create_identity_matrix< xcrsmatrix_t >( a.numRows() );

          if ( timer2_ptr ) {
            space.fence();
            auto duration = timer2_ptr->seconds();
            std::cout
                << "diverg::create_dindex_pairg::create_identity_matrix time: "
                << duration * 1000 << "ms" << std::endl;
          }

          if ( timer2_ptr ) timer2_ptr->reset();

          ai = kokkos_kernels_spadd( a, I );

          if ( timer2_ptr ) {
            space.fence();
            auto duration = timer2_ptr->seconds();
            std::cout
                << "diverg::create_dindex_pairg::kokkos_kernels_spadd time: "
                << duration * 1000 << "ms" << std::endl;
          }
        }  // free: I

        if ( dlo != 0 ) {
          ad = kokkos_kernels_power( a, dlo, timer2_ptr );
        }
      }  // free: a

      aid = kokkos_kernels_power( ai, dup - dlo, timer2_ptr );
    }  // free: ai

    if ( dlo != 0 )
      c = kokkos_kernels_spgemm( ad, aid );
    else
      c = aid;

    if ( timer1_ptr ) {
      space.fence();
      auto duration = timer1_ptr->seconds();
      std::cout << "diverg::create_dindex_pairg time: " << duration * 1000
                << "ms" << std::endl;
    }
  }  // free: ad, aid

  if ( verbose > 1 ) diverg::print( c );

  return c;
}

template< typename TOrdinal, typename TSize, typename TScalar,
          typename TSparseConfig, typename TGraph >
void
benchmark_dindex_graph( TSparseConfig config, TGraph const& graph, int dlo,
                        int dup, typename TGraph::rank_type start_rank,
                        typename TGraph::rank_type end_rank, bool run_kokkos,
                        bool run_rspgemm, bool compare, int verbose )
{
  using ordinal_t = TOrdinal;
  using scalar_t = TScalar;
  using config_type = TSparseConfig;
  using execution_space = typename config_type::execution_space;
  using device_t = typename execution_space::device_type;
  using xcrsmatrix_t = KokkosSparse::CrsMatrix< scalar_t, ordinal_t, device_t >;
  using xcrs_host_mirror = typename xcrsmatrix_t::HostMirror;
  using size_type = std::common_type_t< typename xcrsmatrix_t::size_type, TSize >;
  using range_crsmatrix_t = diverg::CRSMatrix< crs_matrix::RangeDynamic, bool, ordinal_t, size_type >;

  Kokkos::Timer timer;

  std::cout << "Distance index constraints: [d=" << dlo << ", D=" << dup << "]"
            << std::endl;

  std::cout << "Node ranks range: [" << start_rank << ", " << end_rank << ")"
            << std::endl;

  std::cout << "Creating adjacency matrix (on host)..." << std::endl;
  timer.reset();

  auto h_a = diverg::util::adjacency_matrix< xcrs_host_mirror >(
      graph, start_rank, end_rank );

  auto duration = timer.seconds();
  std::cout << "diverg::benchmark_dindex_graph::adjacency_matrix time: "
            << duration * 1000 << "ms" << std::endl;

  std::cout << "Adjacency matrix has order " << h_a.numRows() << "x"
            << h_a.numCols() << " and holds " << h_a.nnz()
            << " non-zero elements" << std::endl;

  xcrs_host_mirror h_c;
  if ( run_kokkos ) {
    std::cout << "Benchmarking PairG..." << std::endl;
    std::cout << "Copy adjacency matrix to device..." << std::endl;
    auto a = copy_xcrs< xcrsmatrix_t >( h_a, &timer );
    if ( verbose > 1 ) diverg::print( a );

    std::shared_ptr< Kokkos::Timer > inner_timer = nullptr;
    if ( verbose > 1 ) inner_timer = std::make_shared< Kokkos::Timer >();
    auto c = create_dindex_pairg( a, dlo, dup, verbose, &timer,
                                  inner_timer.get() );
    // execution_space{}.fence(); // no need to fence here if timer1 is not nullptr
    std::cout << "Distance matrix of order " << c.numRows() << "x"
              << c.numCols() << " holds " << c.nnz() << " non-zero elements"
              << std::endl;

    std::cout << "Copy result from device to host..." << std::endl;
    h_c = copy_xcrs< xcrs_host_mirror >( c, &timer );
  }

  if ( run_rspgemm ) {
    std::cout << "Benchmarking DiVerG..." << std::endl;
    std::cout << "Convert adjacency matrix to range CRS format (on host)..."
              << std::endl;

    range_crsmatrix_t ra( h_a );

    duration = timer.seconds();
    std::cout << "diverg::create_dindex::to_rcrsmatrix(h_a) time (host->host): "
              << duration * 1000 << "ms" << std::endl;

    std::shared_ptr< Kokkos::Timer > inner_timer = nullptr;
    if ( verbose > 1 ) inner_timer = std::make_shared< Kokkos::Timer >();
    auto rc = util::create_distance_index( ra, dlo, dup, config, &timer,
                                           inner_timer.get() );
    // execution_space{}.fence(); // no need to fence here since rc is on host

    // TODO: compare h_c and rc
    // if ( compare && run_kokkos && !is_same( rc, h_c ) ) {
    //   std::cout << "[WARN] Distance matrices are not identical " << std::endl;
    // }

    double comp_rate
        = static_cast< double >( rc.nnz() ) / rc.rowMap( rc.numRows() );
    std::cout << "Distance matrix of order " << rc.numRows() << "x"
              << rc.numCols() << " holds " << rc.nnz()
              << " non-zero elements with compression rates " << comp_rate
              << " (" << rc.rowMap( rc.numRows() ) << ")" << std::endl;

    if ( verbose > 1 ) diverg::print( rc, std::string( "RC" ) );
  }
}

template< typename TOrdinal, typename TSize, typename TScalar,
          typename TSparseConfig >
void
benchmark_dindex_random( TSparseConfig config, TOrdinal n, TSize nnz, int dlo,
                         int dup, unsigned int seed, bool run_kokkos,
                         bool run_rspgemm, bool compare, int verbose )
{
  using ordinal_t = TOrdinal;
  using scalar_t = TScalar;
  using config_type = TSparseConfig;
  using execution_space = typename config_type::execution_space;
  using device_t = typename execution_space::device_type;
  using xcrsmatrix_t = KokkosSparse::CrsMatrix< scalar_t, ordinal_t, device_t >;
  using xcrs_host_mirror = typename xcrsmatrix_t::HostMirror;
  using size_type = std::common_type_t< typename xcrsmatrix_t::size_type, TSize >;
  using range_crsmatrix_t = diverg::CRSMatrix< crs_matrix::RangeDynamic, bool, ordinal_t, size_type >;

  Kokkos::Timer timer;

  std::cout << "Distance index constraints: [d=" << dlo << ", D=" << dup << "]"
            << std::endl;

  std::cout << "Random seed: " << seed << std::endl;

  std::cout << "Creating random binary matrix (on device)..." << std::endl;
  timer.reset();

  // TODO: pass seed
  range_crsmatrix_t ra;
  auto a = create_random_binary_matrix< xcrsmatrix_t >( n, nnz, ra /*, seed*/ );

  auto duration = timer.seconds();
  std::cout << "diverg::benchmark_dindex_random::create_random_binary_matrix time: "
            << duration * 1000 << "ms" << std::endl;

  std::cout << "Random binary matrix has order " << a.numRows() << "x"
            << a.numCols() << " and holds " << a.nnz() << " non-zero elements"
            << std::endl;

  xcrs_host_mirror h_c;
  if ( run_kokkos ) {
    std::cout << "Benchmarking PairG..." << std::endl;
    if ( verbose > 1 ) diverg::print( a );

    std::shared_ptr< Kokkos::Timer > inner_timer = nullptr;
    if ( verbose > 1 ) inner_timer = std::make_shared< Kokkos::Timer >();
    auto c = create_dindex_pairg( a, dlo, dup, verbose, &timer,
                                  inner_timer.get() );
    // execution_space{}.fence(); // no need to fence here if timer1 is not nullptr
    std::cout << "Distance matrix of order " << c.numRows() << "x"
              << c.numCols() << " holds " << c.nnz() << " non-zero elements"
              << std::endl;

    std::cout << "Copy result from device to host..." << std::endl;
    h_c = copy_xcrs< xcrs_host_mirror >( c, &timer );
  }

  if ( run_rspgemm ) {
    std::cout << "Benchmarking DiVerG..." << std::endl;

    std::shared_ptr< Kokkos::Timer > inner_timer = nullptr;
    if ( verbose > 1 ) inner_timer = std::make_shared< Kokkos::Timer >();
    auto rc = util::create_distance_index( ra, dlo, dup, config, &timer,
                                           inner_timer.get() );
    // execution_space{}.fence(); // no need to fence here since rc is on host

    // TODO: compare h_c and rc
    // if ( compare && run_kokkos && !is_same( rc, h_c ) ) {
    //   std::cout << "[WARN] Distance matrices are not identical " << std::endl;
    // }

    double comp_rate
        = static_cast< double >( rc.nnz() ) / rc.rowMap( rc.numRows() );
    std::cout << "Distance matrix of order " << rc.numRows() << "x"
              << rc.numCols() << " holds " << rc.nnz()
              << " non-zero elements with compression rates " << comp_rate
              << " (" << rc.rowMap( rc.numRows() ) << ")" << std::endl;

    if ( verbose > 1 ) diverg::print( rc, std::string( "RC" ) );
  }
}

template< typename TOrdinal, typename TSize, typename TScalar, typename TAcc,
          typename TExecSpace, typename TPartition, typename TGrid >
void
benchmark_dindex( Options< TOrdinal, TSize > opts,
                  SparseConfig< TGrid, TAcc, TPartition, TExecSpace > config )
{
  if ( opts.graph_path.empty() ) {
    benchmark_dindex_random< TOrdinal, TSize, TScalar >(
        config, opts.n, opts.nnz, opts.dlo, opts.dup, opts.seed, opts.run_kokkos,
        opts.run_rspgemm, opts.compare, opts.verbose );
  }
  else {
    using graph_type = gum::SeqGraph< gum::Succinct >;

    auto graph = load_graph< graph_type >( opts.graph_path, opts.format );
    auto rank_range = region_nodes_rank_range( graph, opts.reg_name, opts.seg_name );
    benchmark_dindex_graph< TOrdinal, TSize, TScalar >(
        config, graph, opts.dlo, opts.dup, rank_range.first, rank_range.second,
        opts.run_kokkos, opts.run_rspgemm, opts.compare, opts.verbose );
  }
}

template< typename TOrdinal, typename TSize, typename TScalar, typename TAcc,
          typename TExecSpace, typename TPartition >
void
benchmark_dindex( Options< TOrdinal, TSize > opts, TExecSpace space,
                  TPartition partition )
{
  if ( opts.grid == "auto" ) {
    std::cout << "Grid: Auto" << std::endl;
    using grid_type = grid::Auto;
    using config_type
        = SparseConfig< grid_type, TAcc, TPartition, TExecSpace >;
    benchmark_dindex< TOrdinal, TSize, TScalar >( opts, config_type{} );
  }
  else if ( opts.grid == "suggested" ) {
    std::cout << "Grid: Suggested" << std::endl;
    using grid_type = grid::Suggested;
    using config_type
        = SparseConfig< grid_type, TAcc, TPartition, TExecSpace >;
    benchmark_dindex< TOrdinal, TSize, TScalar >( opts, config_type{} );
  }
  else {
    std::cout << "Grid: RunTime" << std::endl;
    using grid_type = grid::RunTime;
    using config_type
        = SparseConfig< grid_type, TAcc, TPartition, TExecSpace >;
    config_type config;
    config.grid.m_team_work_size = RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE;
    try {
      std::cout << "Parsing provided grid dimension..." << std::endl;
      std::istringstream iss( opts.grid );
      char delimiter;
      if ( !( iss >> config.grid.m_team_size ) )
        throw std::invalid_argument( "parsing team size" );
      if ( !( iss >> delimiter ) )
        throw std::invalid_argument( "parsing first delimiter" );
      if ( !( iss >> config.grid.m_vector_size ) )
        throw std::invalid_argument( "parsing vector size" );
      if ( iss >> delimiter )
        iss >> config.grid.m_team_work_size;
    }
    catch ( std::invalid_argument& e ) {
      std::cerr << "[ERROR] Invalid grid option '" << opts.grid
                << "': " << e.what() << std::endl;
      exit( EXIT_FAILURE );
    }
    if ( config.grid.team_size() == 0 || config.grid.vector_size() == 0
         || config.grid.team_work_size() == 0 ) {
      std::cerr << "[ERROR] Grid elements should be non-zero" << std::endl;
      exit( EXIT_FAILURE );
    }
    std::cout << "Grid< team size, vector size, team work size >: <"
              << config.grid.team_size() << ", " << config.grid.vector_size()
              << ", " << config.grid.team_work_size() << ">" << std::endl;
    benchmark_dindex< TOrdinal, TSize, TScalar >( opts, config );
  }
}

template< typename TOrdinal, typename TSize, typename TScalar, typename TAcc,
          typename TExecSpace >
void
benchmark_dindex( Options< TOrdinal, TSize > opts, TExecSpace space )
{
  if ( opts.partition == "thread-sequential" ) {
    std::cout << "Partition: Thread-Sequential" << std::endl;
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >(
        opts, space, ThreadSequentialPartition{} );
  }
  else if ( opts.partition == "team-sequential" ) {
    std::cout << "Partition: Team-Sequential" << std::endl;
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >(
        opts, space, TeamSequentialPartition{} );
  }
  else if ( opts.partition == "thread-parallel" ) {
    std::cout << "Partition: Thread-Parallel" << std::endl;
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >(
        opts, space, ThreadParallelPartition{} );
  }
  else if ( opts.partition == "thread-range-parallel" ) {
    std::cout << "Partition: Thread-Range-Parallel" << std::endl;
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >(
        opts, space, ThreadRangeParallelPartition{} );
  }
  else
    throw std::runtime_error( "[ERROR] Unknown partition scheme '"
                              + opts.partition + "'" );
}

template< typename TOrdinal, typename TSize, typename TScalar, typename TAcc >
void
benchmark_dindex( Options< TOrdinal, TSize > opts )
{
  if ( opts.host ) {
    using execution_space = Kokkos::DefaultHostExecutionSpace;
    execution_space space;
    std::cout << "Execution space: Host" << std::endl;
    std::cout << "Execution space concurrency: "
              << execution_space::concurrency() << std::endl;
    space.print_configuration( std::cout, opts.verbose );
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >( opts, space );
  }
  else {
    using execution_space = Kokkos::DefaultExecutionSpace;
    execution_space space;
#if defined(KOKKOS_ENABLE_CUDA)
    if constexpr ( std::is_same< execution_space, Kokkos::Cuda >::value ) {
      std::cout << "Execution space: Cuda" << std::endl;
    }
    else {
      std::cout
          << "Execution space: Default (Kokkos::DefaultExecutionSpace != Cuda)"
          << std::endl;
    }
#else
    std::cout << "Execution space: Default" << std::endl;
#endif
    std::cout << "Execution space concurrency: "
              << execution_space::concurrency() << std::endl;
    space.print_configuration( std::cout, opts.verbose );
    benchmark_dindex< TOrdinal, TSize, TScalar, TAcc >( opts, space );
  }
}

template< typename TOrdinal, typename TSize, typename TScalar >
void
benchmark_dindex( Options< TOrdinal, TSize > opts )
{
  if ( opts.l1size == 1024 ) {
    std::cout << "L1 size: 1024" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 1024 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else if ( opts.l1size == 2048 ) {
    std::cout << "L1 size: 2048" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 2048 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else if ( opts.l1size == 4096 ) {
    std::cout << "L1 size: 4096" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 4096 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else if ( opts.l1size == 8192 ) {
    std::cout << "L1 size: 8192" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 8192 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else if ( opts.l1size == 16384 ) {
    std::cout << "L1 size: 16384" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 16384 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else if ( opts.l1size == 32768 ) {
    std::cout << "L1 size: 32768" << std::endl;
    using accumulator_type = HBitVectorAccumulator< 32768 >;
    benchmark_dindex< TOrdinal, TSize, TScalar, accumulator_type >( opts );
  }
  else
    throw std::runtime_error( "[ERROR] Invalid L1 size '"
                              + std::to_string( opts.l1size ) + "'" );
}

template< typename TOrdinal = int32_t, typename TSize = uint64_t >
struct Options {
  /* === TYPE MEMBERS === */
  using ordinal_type = TOrdinal;
  using size_type = TSize;
  /* === DATA MEMBERS === */
  std::string graph_path;
  std::string format;
  TOrdinal n;
  TSize nnz;
  unsigned int seed;
  unsigned int dlo;
  unsigned int dup;
  std::string seg_name;
  std::string reg_name;
  std::string partition;
  std::string grid;
  unsigned int l1size;
  bool run_kokkos;
  bool run_rspgemm;
  bool host;
  bool compare;
  int verbose;
  /* === LIFECYCLE === */
  Options()
      : graph_path( "" ), format( "" ), n( 0 ), nnz( 0 ), dlo( 0 ), dup( 0 ),
        seg_name( "" ), reg_name( "" ), partition( "team-sequential" ),
        grid( "auto" ), l1size( 8192 ), run_kokkos( true ),
        run_rspgemm( true ), host( false ), compare( false ), verbose( 0 )
  { }
};

template< typename TOrdinal = Options<>::ordinal_type,
          typename TSize = Options<>::size_type >
Options< TOrdinal, TSize >
parse_arguments( int argc, char* argv[] )
{
  Options< TOrdinal, TSize > opts;

  // Parse command line arguments
  for ( int i = 1; i < argc; i++ ) {
    if ( ( strcmp( argv[ i ], "-g" ) == 0 ) || ( strcmp( argv[ i ], "--graph" ) == 0 ) ) {
      opts.graph_path = argv[ ++i ];
      std::cout << "Parameter: graph <- " << opts.graph_path << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-f" ) == 0 ) || ( strcmp( argv[ i ], "--format" ) == 0 ) ) {
      opts.format = argv[ ++i ];
      std::cout << "Parameter: format <- " << opts.format << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-n" ) == 0 ) || ( strcmp( argv[ i ], "--order" ) == 0 ) ) {
      opts.n = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: n <- " << opts.n << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-z" ) == 0 ) || ( strcmp( argv[ i ], "--nnz" ) == 0 ) ) {
      opts.nnz = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: nnz <- " << opts.nnz << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-S" ) == 0 ) || ( strcmp( argv[ i ], "--seed" ) == 0 ) ) {
      opts.seed = atoi( argv[ ++i ] );
      std::cout << "Parameter: seed <- " << opts.seed << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-d" ) == 0 ) || ( strcmp( argv[ i ], "--dist-lower" ) == 0 ) ) {
      opts.dlo = atoi( argv[ ++i ] );
      std::cout << "Parameter: d <- " << opts.dlo << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-D" ) == 0 ) || ( strcmp( argv[ i ], "--dist-upper" ) == 0 ) ) {
      opts.dup = atoi( argv[ ++i ] );
      std::cout << "Parameter: D <- " << opts.dup << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-K" ) == 0 ) || ( strcmp( argv[ i ], "--disable-kokkos" ) == 0 ) ) {
      opts.run_kokkos = false;
      std::cout << "Parameter: disable-kokkos <- "
                << ( opts.run_kokkos ? "false" : "true" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-R" ) == 0 ) || ( strcmp( argv[ i ], "--disable-rspgemm" ) == 0 ) ) {
      opts.run_rspgemm = false;
      std::cout << "Parameter: disable-rspgemm <- "
                << ( opts.run_rspgemm ? "false" : "true" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-s" ) == 0 ) || ( strcmp( argv[ i ], "--segment-name" ) == 0 ) ) {
      opts.seg_name = argv[ ++i ];
      std::cout << "Parameter: segment-name <- " << opts.seg_name << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-r" ) == 0 ) || ( strcmp( argv[ i ], "--region-name" ) == 0 ) ) {
      opts.reg_name = argv[ ++i ];
      std::cout << "Parameter: region-name <- " << opts.reg_name << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-p" ) == 0 ) || ( strcmp( argv[ i ], "--partition" ) == 0 ) ) {
      opts.partition = argv[ ++i ];
      std::cout << "Parameter: partition <- " << opts.partition << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-G" ) == 0 ) || ( strcmp( argv[ i ], "--grid" ) == 0 ) ) {
      opts.grid = argv[ ++i ];
      std::cout << "Parameter: grid <- " << opts.grid << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-l" ) == 0 ) || ( strcmp( argv[ i ], "--l1-size" ) == 0 ) ) {
      opts.l1size = pow( 2, atoi( argv[ ++i ] ) );
      std::cout << "Parameter: l1size <- " << opts.l1size << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-H" ) == 0 ) || ( strcmp( argv[ i ], "--host" ) == 0 ) ) {
      opts.host = true;
      std::cout << "Parameter: host <- "
                << ( opts.host ? "true" : "false" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-c" ) == 0 ) || ( strcmp( argv[ i ], "--compare" ) == 0 ) ) {
      opts.compare = true;
      std::cout << "Parameter: compare <- "
                << ( opts.compare ? "true" : "false" ) << std::endl;
    }
    else if ( ( strcmp( argv[ i ], "-v" ) == 0 ) || ( strcmp( argv[ i ], "--verbose" ) == 0 ) ) {
      opts.verbose++;
    }
    else if ( ( strcmp( argv[ i ], "-h" ) == 0 ) || ( strcmp( argv[ i ], "--help" ) == 0 ) ) {
      std::cout
          << "Options:\n"
          << "  * Adjacency matrix as input (overrides -n and -z)\n"
          << "    --graph (-g) <path>:       graph file path\n"
          << "    --format (-f) <str>:       specify graph file format when it is not\n"
          << "                               automatically detected; e.g. when GFA has no header\n"
          << "                               (formats: 'gfa', 'gfa1', 'gfa2')\n"
          << "    --region-name (-r) <str>:  specify target region by path id (embedded in graph)\n"
          << "                               (overrides -s)\n"
          << "    --segment-name (-s) <str>: specify target region by the name of its first node\n"
          << "                               or segment (default: first node/segment in the graph)\n"
          << "  * Random matrix as input\n"
          << "    --order (-n) <int>:        exponent num, determines number of rows 2^num\n"
          << "                               (default: 2^12 = 4096 or 2^((z+1)/2) if -z given)\n"
          << "    --nnz (-z) <int>:          exponent num, determines total matrix size 2^num\n"
          << "                               (default: 2^22 = 4096*1024 or 2^n*1024 if -n provided)\n"
          << "    --rng-seed (-S) <int>:     Seed for random number generation (default: random)\n"
          << "  --dist-lower (-d) <int>:     lower bound of outer distance, i.e. fragment size\n"
          << "                               (default: 0)\n"
          << "  --dist-upper (-D) <int>:     upper bound of outer distance, i.e. fragment size\n"
          << "                               (required)\n"
          << "  --partition (-p) <str>:      Partitioning schemes: 'team-sequential' (default),\n"
          << "                               'thread-sequential', 'thread-parallel', and\n"
          << "                               'thread-range-paralell'\n"
          << "  --grid (-G) <str>            Execution grid options: 'auto' (auto-tuning; default),\n"
          << "                               'suggested' (based on row density of input)\n"
          << "  --grid (-G) <int,int[,int]>  '<team size>,<vector size>[,<team work size>]'\n"
          << "                               in which default team work size is "
          << RCRS_BENCHMARK_DEFAULT_TEAM_WORK_SIZE_STR << "\n"
          << "                               (any non-numeric, non-whitespace delimiter works)\n"
          << "                               examples: '-G 32,4,16' or '-G 16x16'\n"
          << "  --l1-size (-l) <int>:        exponent num, determines L1 size of accumulator 2^num;\n"
          << "                               with constraints: 10<=num<=15 (default: 2^13 = 1KB)\n"
          << "  --host (-H):                 run on host (default: false)\n"
          << "  --disable-kokkos (-K):       disable benchmarking with Kokkos SpGEMM (PairG)\n"
          << "                               (default: false)\n"
          << "  --disable-rspgemm (-R):      disable benchmarking with rSpGEMM (DiVerG)\n"
          << "                               (default: false)\n"
          << "  --compare (-c):              Compare resulting indices build by DiVerG and PairG\n"
          << "                               given that -K is not set (default: false)\n"
          << "  --verbose (-v):              increase verbose level by one (default: 0)\n"
          << "  --help (-h):                 print this message"
          << std::endl;
      exit( EXIT_SUCCESS );
    }
    else {
      std::cerr << "[ERROR] Unknown argument '" << argv[ i ] << "'" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( opts.verbose ) std::cout << "Verbose level: " << opts.verbose << std::endl;

  return opts;
}

template< typename TOrdinal, typename TSize >
void
check_options( Options< TOrdinal, TSize >& opts )
{
  if ( opts.graph_path.empty() ) {
    if ( opts.nnz == 0 ) {
      if ( opts.n == 0 )  // If both opts.n and opts.nnz are not provided:
        opts.nnz = pow( 2, 22 );
      else  // If opts.n is provided but opts.nnz is not:
        opts.nnz = opts.n * 1024;
      std::cout << "[DEFAULT] nnz = " << opts.nnz << std::endl;
    }

    // If opts.n is not provided, set order to 2^((log2(opts.nnz)+1)/2).
    if ( opts.n == 0 ) {
      int power = 1;
      auto temp = opts.nnz;
      while ( ( temp /= 2 ) ) ++power;
      opts.n = pow( 2, power / 2 );
      std::cout << "[DEFAULT] n = " << opts.n << std::endl;
    }

    if ( opts.n < 0 /* || opts.nnz < 0 */ ) {
      std::cerr << "[ERROR] nnz must be greater than zero" << std::endl;
      exit( EXIT_FAILURE );
    }

    if ( opts.nnz / opts.n > static_cast< TSize >( opts.n ) ) {
      std::cerr << "[ERROR] Non-zero values cannot be fit (nnz > n*n)"
                << std::endl;
      exit( EXIT_FAILURE );
    }
  }
  else {
    if ( !gum::util::readable( opts.graph_path ) ) {
      std::cerr << "[ERROR] Input file is missing or not readable"
                << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  if ( opts.dup == 0 ) {
    std::cerr << "[ERROR] Too few arguments (required argument: '-D')"
              << std::endl;
    exit( EXIT_FAILURE );
  }

  if ( opts.dup < opts.dlo ) {
    std::cerr << "[WARN] Lower bound of distance constraints is larger than "
                 "upper bound"
              << std::endl;
    std::swap( opts.dlo, opts.dup );
    std::cerr << "[WARN] Swapped distance constraints bounds: (" << opts.dlo
              << ", " << opts.dup << ")" << std::endl;
  }

  if ( !opts.run_kokkos && !opts.run_rspgemm ) {
    std::cerr << "[ERROR] Nothing to do!" << std::endl;
    exit( EXIT_FAILURE );
  }
}

int
main( int argc, char* argv[] )
{
  using ordinal_type = int32_t;
  using size_type = uint64_t;
  using scalar_type = int;

  auto opts = parse_arguments< ordinal_type, size_type >( argc, argv );
  check_options( opts );

  Kokkos::initialize( argc, argv );

  benchmark_dindex< ordinal_type, size_type, scalar_type >( opts );

  Kokkos::finalize();

  return EXIT_SUCCESS;
}
