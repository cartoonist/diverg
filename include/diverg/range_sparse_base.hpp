/**
 *    @file  range_sparse_base.hpp
 *   @brief  Base header for range sparse matrix operations module
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Thu Sep 07, 2023  21:37
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_RANGE_SPARSE_BASE_HPP_
#define DIVERG_RANGE_SPARSE_BASE_HPP_

#include <cstddef>
#include <stdexcept>
#include <algorithm>

#include <Kokkos_Core.hpp>
#include <Kokkos_NestedSort.hpp>


namespace diverg {
  // Accumulator tag
  template< typename TSpec >
  struct Accumulator {
    using spec_type = TSpec;
  };

  // Accumulator specialisation tag
  struct BTreeTag {};

  template< unsigned int TL1Size >
  struct HBitVectorTag {
    static constexpr const unsigned int L1SizeParam = TL1Size;
  };

  using NoAccumulator = Accumulator< void >;
  using BTreeAccumulator = Accumulator< BTreeTag >;

  template< unsigned int TL1Size=2048 >
  using HBitVectorAccumulator = Accumulator< HBitVectorTag< TL1Size > >;

  // Partition tag
  template< typename TSpec >
  struct ExecPartition {
    using type = TSpec;
  };

  // Partition specialisation tag
  struct ThreadRangePolicyTag {};
  struct ThreadSequentialTag {};
  struct TeamSequentialTag {};
  struct ThreadParallelTag {};
  struct ThreadRangeParallelTag {};
  struct TeamFlatParallelTag {};

  using ThreadRangePolicyPartition = ExecPartition< ThreadRangePolicyTag >;
  using ThreadSequentialPartition = ExecPartition< ThreadSequentialTag >;
  using TeamSequentialPartition = ExecPartition< TeamSequentialTag >;
  using ThreadParallelPartition = ExecPartition< ThreadParallelTag >;
  using ThreadRangeParallelPartition = ExecPartition< ThreadRangeParallelTag >;
  using TeamFlatParallelPartition = ExecPartition< TeamFlatParallelTag >;

  // Supported execution space by accumulators
  template< typename TAccumulator >
  struct AccumulatorExecSpace {
    using type = Kokkos::DefaultExecutionSpace;      // By default run on device
  };

  template<>
  struct AccumulatorExecSpace< BTreeAccumulator > {
    using type = Kokkos::DefaultHostExecutionSpace;  // Can only be executed on host
  };

  // Default partitioning by accumulators
  template< typename TAccumulator >
  struct AccumulatorDefaultPartition;

  template< >
  struct AccumulatorDefaultPartition< BTreeAccumulator >
  {
    using type = ThreadRangePolicyPartition;
  };

  template< unsigned int TL1Size >
  struct AccumulatorDefaultPartition< HBitVectorAccumulator< TL1Size > >
  {
    using type = TeamSequentialPartition;
  };

  namespace grid {
    template< int TTeamSize, int TVectorSize, int TTeamWorkSize=0 >
    struct Fixed {};
    struct Auto {};
    struct Suggested {};
    struct RunTime {};
  }

  struct _ExecGridBase {
    static constexpr inline int
    row_density( const std::size_t nnz, const std::size_t nr )
    {
      if ( nr > 0 ) return nnz / double( nr ) + 0.5;
      return 1;
    }
  };

  template< typename TExecSpace, typename TSpec >
  struct ExecGrid;

  template< typename TExecSpace,
            int TTeamSize,
            int TVectorSize,
            int TTeamWorkSize >
  struct ExecGrid< TExecSpace, grid::Fixed< TTeamSize, TVectorSize, TTeamWorkSize > >
      : public _ExecGridBase {
    static constexpr inline int
    vector_size( const int=0 )
    {
      return TVectorSize;
    }

    static constexpr inline int
    vector_size( const std::size_t, const std::size_t )
    {
      return TVectorSize;
    }

    static constexpr inline int
    team_size( const int=0 )
    {
      return TTeamSize;
    }

    static constexpr inline int
    team_size( const std::size_t, const std::size_t )
    {
      return TTeamSize;
    }

    static constexpr inline int
    team_work_size( const int=0 )
    {
      return TTeamWorkSize;
    }

    static constexpr inline int
    team_work_size( const std::size_t, const std::size_t )
    {
      return TTeamWorkSize;
    }
  };

  /**
   *   @brief Suggested grid dimensions as a form of config class.
   */
  template< typename TExecSpace >
  struct ExecGrid< TExecSpace, grid::Suggested >
      : public _ExecGridBase {
    static constexpr inline int
    vector_size( const int=0 )
    {
      return 1;
    }

    static constexpr inline int
    vector_size( const std::size_t, const std::size_t )
    {
      return 1;
    }

    static constexpr inline int
    team_size( const int=0 )
    {
      return 1;
    }

    static constexpr inline int
    team_size( const std::size_t, const std::size_t )
    {
      return 1;
    }

    static constexpr inline int
    team_work_size( const int=0 )
    {
      return 16;
    }

    static constexpr inline int
    team_work_size( const std::size_t, const std::size_t )
    {
      return 16;
    }
  };

  template< typename TExecSpace >
  struct ExecGrid< TExecSpace, grid::Auto >
      : public _ExecGridBase {
    static constexpr inline auto
    vector_size( const int=0 )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const std::size_t, const std::size_t )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const int=0 )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const std::size_t, const std::size_t )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline int
    team_work_size( const int=0 )
    {
      return 16;
    }

    static constexpr inline int
    team_work_size( const std::size_t, const std::size_t )
    {
      return 16;
    }
  };

  template< typename TExecSpace >
  struct ExecGrid< TExecSpace, grid::RunTime >
      : public _ExecGridBase {
    int m_vector_size;
    int m_team_size;
    int m_team_work_size;

    ExecGrid( const int vs=0, const int ts=0, const int tws=0 )
        : m_vector_size( vs ), m_team_size( ts ), m_team_work_size( tws )
    { }

    inline int
    vector_size( const int=0 )
    {
      return this->m_vector_size;
    }

    inline int
    vector_size( const std::size_t, const std::size_t )
    {
      return this->m_vector_size;
    }

    inline int
    team_size( const int=0 )
    {
      return this->m_team_size;
    }

    inline int
    team_size( const std::size_t, const std::size_t )
    {
      return this->m_team_size;
    }

    inline int
    team_work_size( const int=0 )
    {
      return this->m_team_work_size;
    }

    inline int
    team_work_size( const std::size_t, const std::size_t )
    {
      return this->m_team_work_size;
    }
  };

  #if defined(KOKKOS_ENABLE_CUDA)
  template< >
  struct ExecGrid< Kokkos::Cuda, grid::Suggested >
      : public _ExecGridBase {
    /* === MEMBER TYPES === */
    using base_type = _ExecGridBase;
    /* === STATIC MEMBERS === */
    static constexpr const int MAX_VECTOR_SIZE = 32;
    /* === STATIC METHODS === */
    static constexpr inline int
    vector_size( const int rdense )
    {
      int vsize = rdense;
      if ( vsize < 3 ) {
        vsize = 2;
      } else if ( vsize <= 6 ) {
        vsize = 4;
      } else if ( vsize <= 12 ) {
        vsize = 8;
      } else if ( vsize <= 24 ) {
        vsize = 16;
      } else if ( vsize <= 48 ) {
        vsize = 32;
      } else {
        vsize = 64;
      }
      vsize = Kokkos::min( vsize, MAX_VECTOR_SIZE );
      return vsize;
    }

    static constexpr inline int
    vector_size( const std::size_t nnz, const std::size_t nr )
    {
      return ExecGrid::vector_size( base_type::row_density( nnz, nr ) );
    }

    static constexpr inline int
    team_size( const int rdense )
    {
      // TODO: where this is used, tune the target value for
      // threads per block (but 256 is probably OK for CUDA and HIP)
      return 256 / ExecGrid::vector_size( rdense );
    }

    static constexpr inline int
    team_size( const std::size_t nnz, const std::size_t nr )
    {
      return ExecGrid::team_size( base_type::row_density( nnz, nr ) );
    }

    static constexpr inline int
    team_work_size( const int rdense )
    {
      return ExecGrid::team_size( rdense );
    }

    static constexpr inline int
    team_work_size( const std::size_t nnz, const std::size_t nr )
    {
      return ExecGrid::team_work_size( base_type::row_density( nnz, nr ) );
    }
  };

  template< >
  struct ExecGrid< Kokkos::Cuda, grid::Auto >
      : public _ExecGridBase {
    /* === MEMBER TYPES === */
    using base_type = _ExecGridBase;
    /* === STATIC MEMBERS === */
    static constexpr const int MAX_VECTOR_SIZE = 32;
    /* === STATIC METHODS === */
    static constexpr inline auto
    vector_size( const int=0 )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    vector_size( const std::size_t, const std::size_t )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const int=0 )
    {
      return Kokkos::AUTO;
    }

    static constexpr inline auto
    team_size( const std::size_t, const std::size_t )
    {
      return Kokkos::AUTO;
    }

    static inline int
    team_work_size( const int rdense )
    {
      int vsize = rdense;
      if ( vsize < 3 ) {
        vsize = 2;
      } else if ( vsize <= 6 ) {
        vsize = 4;
      } else if ( vsize <= 12 ) {
        vsize = 8;
      } else if ( vsize <= 24 ) {
        vsize = 16;
      } else if ( vsize <= 48 ) {
        vsize = 32;
      } else {
        vsize = 64;
      }
      vsize = Kokkos::min( vsize, MAX_VECTOR_SIZE );
      return 256 / vsize;
    }

    static inline int
    team_work_size( const std::size_t nnz, const std::size_t nr )
    {
      return ExecGrid::team_work_size( base_type::row_density( nnz, nr ) );
    }
  };
  #endif

  /* === ExecGrid meta-functions === */

  template< typename TExecSpace, typename TSpec >
  struct GetExecGrid {
    using type = ExecGrid< TExecSpace, TSpec >;
  };

  template< typename TExecSpace, typename TSpec >
  using ExecGridType = typename GetExecGrid< TExecSpace, TSpec >::type;

  template< typename TTargetExecSpace,
            typename TExecSpace1, typename TGridSpec1,
            typename TExecSpace2=void, typename TGridSpec2=grid::Auto >
  struct MatchingGridSpec {
    using type = void;
  };

  /* NOTE: Keep default argument types in sync with `MatchingGridSpec` */
  template< typename TTargetExecSpace,
            typename TExecSpace1, typename TGridSpec1,
            typename TExecSpace2=void, typename TGridSpec2=grid::Auto >
  using MatchingGridSpecType =
      typename MatchingGridSpec< TTargetExecSpace, TExecSpace1, TGridSpec1,
                                 TExecSpace2, TGridSpec2 >::type;

  template< typename TExecSpace, typename TExecSpace2,
            typename TGridSpec1, typename TGridSpec2 >
  struct MatchingGridSpec< TExecSpace, TExecSpace, TGridSpec1,
                           TExecSpace2, TGridSpec2 >
  {
    using type = TGridSpec1;
  };

  template< typename TExecSpace, typename TExecSpace1,
            typename TGridSpec1, typename TGridSpec2 >
  struct MatchingGridSpec< TExecSpace, TExecSpace1, TGridSpec1,
                           TExecSpace, TGridSpec2 >
  {
    using type = TGridSpec2;
  };

  template< typename TExecSpace, typename TGridSpec1, typename TGridSpec2 >
  struct MatchingGridSpec< TExecSpace, TExecSpace, TGridSpec1,
                           TExecSpace, TGridSpec2 >
  {
    using type = TGridSpec1;  // prefer the first match
  };

  template< typename TExecSpace, typename TGridSpec1, typename TGridSpec2 >
  struct MatchingGridSpec< TExecSpace, TExecSpace, TGridSpec1,
                           void, TGridSpec2 > {
    using type = TGridSpec1;
  };

  template< typename TExecSpace, typename TExecSpace1,
            typename TGridSpec1, typename TGridSpec2 >
  struct MatchingGridSpec< TExecSpace, TExecSpace1, TGridSpec1,
                           void, TGridSpec2 > {
    using type = TGridSpec2;  // use default
  };

  /* === End of ExecGrid meta-functions === */

  // Configuration tag
  template< typename TGridSpec,
            typename TAccumulator,
            typename TPartition =
                typename AccumulatorDefaultPartition< TAccumulator >::type,
            typename TExecSpace = typename AccumulatorExecSpace< TAccumulator >::type >
  struct SparseConfig {
    using partition_type = TPartition;
    using accumulator_type = TAccumulator;
    using execution_space = TExecSpace;
    using grid_type = ExecGridType< execution_space, TGridSpec >;
    /* === MEMBERS === */
    partition_type   part;
    accumulator_type accm;
    execution_space  space;
    grid_type        grid;
  };

  using DefaultSparseConfiguration = SparseConfig< grid::Auto, HBitVectorAccumulator<> >;

  template< typename TRCRSMatrix, typename TExecSpace >
  struct SparseRangeHandle {
    /* === MEMBERS TYPES === */
    using ordinal_type = typename TRCRSMatrix::ordinal_type;
    using size_type = typename TRCRSMatrix::size_type;
    using execution_space = TExecSpace;
    using entries_device_view_type =
        typename TRCRSMatrix::template entries_device_view_type< execution_space >;
    /* === LIFECYCLE === */
    SparseRangeHandle() = delete;

    SparseRangeHandle( TRCRSMatrix const& a, TRCRSMatrix const& b, TExecSpace={} )
        : a_ncols( a.numCols() ), b_ncols( b.numCols() ), a_nnz( a.nnz() ),
          b_nnz( b.nnz() ), c_bandwidth( 0 ), c_min_col_index()
    { }

    SparseRangeHandle( ordinal_type ncols_a, size_type nnz_a,
                       ordinal_type ncols_b, size_type nnz_b, TExecSpace={} )
        : a_ncols( ncols_a ), b_ncols( ncols_b ), a_nnz( nnz_a ),
          b_nnz( nnz_b ), c_bandwidth( 0 ), c_min_col_index()
    { }
    /* === METHODS === */
    inline void init_c_min_col_index( size_type nrows )
    {
      this->c_min_col_index = entries_device_view_type( "c_min_col_index", nrows );
    }
    /* === DATA MEMBERS === */
    ordinal_type a_ncols;
    ordinal_type b_ncols;
    size_type a_nnz;
    size_type b_nnz;
    size_type c_bandwidth;
    entries_device_view_type c_min_col_index;
  };

  template<
      typename TRowMapView,
      typename TEntriesView,
      typename TExecSpace=typename TRowMapView::execution_space >
  struct SortEntriesFunctor {
    /* === TYPE MEMBERS === */
    using policy_type = Kokkos::RangePolicy< TExecSpace >;
    using size_type = typename TRowMapView::non_const_value_type;
    /* === STATIC ASSERTS === */
    static_assert(
        std::is_same< typename TRowMapView::execution_space,
                      TExecSpace >::value,
        "execution space parameter is incompatible with input views" );
    static_assert(
        std::is_same< typename TEntriesView::execution_space,
                      TExecSpace >::value,
        "execution space parameter is incompatible with input views" );
    /* === DATA MEMBERS === */
    TRowMapView row_map;
    TEntriesView entries;
    /* === LIFE CYCLE === */
    SortEntriesFunctor( TRowMapView r, TEntriesView e )
      : row_map( r ), entries( e )
    { }
    /* === METHODS === */
    inline auto
    policy( const size_type nrows ) const
    {
      return policy_type( 0, nrows );
    }
    /* === OPERATORS === */
    inline void
    operator()( const uint64_t i ) const
    {
      auto begin = this->entries.data() + this->row_map( i );
      auto end = this->entries.data() + this->row_map( i + 1 );
      std::sort( begin, end );
    }
  };

#if defined(KOKKOS_ENABLE_CUDA)
  template< typename TRowMapView, typename TEntriesView >
  struct SortEntriesFunctor< TRowMapView, TEntriesView, Kokkos::Cuda > {
    /* === TYPE MEMBERS === */
    using execution_space = Kokkos::Cuda;
    using policy_type = Kokkos::TeamPolicy< execution_space >;
    using member_type = typename policy_type::member_type;
    using size_type = typename TRowMapView::non_const_value_type;
    /* === STATIC ASSERTS === */
    static_assert(
        std::is_same< typename TRowMapView::execution_space,
                      execution_space >::value,
        "execution space parameter is incompatible with input views" );
    static_assert(
        std::is_same< typename TEntriesView::execution_space,
                      execution_space >::value,
        "execution space parameter is incompatible with input views" );
    /* === DATA MEMBERS === */
    TRowMapView row_map;
    TEntriesView entries;
    /* === LIFE CYCLE === */
    SortEntriesFunctor( TRowMapView r, TEntriesView e )
      : row_map( r ), entries( e )
    { }
    /* === METHODS === */
    inline auto
    policy( const size_type nrows ) const
    {
      return policy_type( nrows, Kokkos::AUTO );
    }
    /* === OPERATORS === */
    KOKKOS_INLINE_FUNCTION void
    operator()( member_type const& tm ) const
    {
      auto i = tm.league_rank();
      auto l = this->row_map( i );
      auto u = this->row_map( i + 1 );
      assert( u >= l );
      if ( u != l ) {
        auto subview = Kokkos::subview( this->entries, Kokkos::pair( l, u ) );
        Kokkos::Experimental::sort_team( tm, subview );
      }
    }
  };
#endif

  /* borrowed from kokkos-kernel */
  template< typename TExecSpace >
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             TExecSpace )
  {
    throw std::runtime_error( "`get_suggested_vector_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int
  get_suggested_vector_size( const std::size_t nr,
                             const std::size_t nnz,
                             Kokkos::Cuda )
  {
    int suggested_vector_size_ = 1;
    int max_vector_size        = 32;
    if ( nr > 0 ) suggested_vector_size_ = nnz / double( nr ) + 0.5;
    if ( suggested_vector_size_ < 3 ) {
      suggested_vector_size_ = 2;
    } else if ( suggested_vector_size_ <= 6 ) {
      suggested_vector_size_ = 4;
    } else if ( suggested_vector_size_ <= 12 ) {
      suggested_vector_size_ = 8;
    } else if ( suggested_vector_size_ <= 24 ) {
      suggested_vector_size_ = 16;
    } else if ( suggested_vector_size_ <= 48 ) {
      suggested_vector_size_ = 32;
    } else {
      suggested_vector_size_ = 64;
    }
    if ( suggested_vector_size_ > max_vector_size )
      suggested_vector_size_ = max_vector_size;
    return suggested_vector_size_;
  }
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::OpenMP )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::Serial )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int
  get_suggested_vector_size( const std::size_t,
                             const std::size_t,
                             Kokkos::Threads )
  {
    //int suggested_vector_size_ = 1;
    //int max_vector_size        = 1;
    return /*suggested_vector_size_*/ 1;
  }
  #endif

  template< typename TExecSpace, int TVectorSize >
  struct SuggestedTeamSize;

  template< typename TExecSpace >
  inline int get_suggested_team_size( const int vector_size, TExecSpace )
  {
    throw std::runtime_error( "`get_suggested_team_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Cuda )
  {
    // TODO: where this is used, tune the target value for
    // threads per block (but 256 is probably OK for CUDA and HIP)
    return 256 / vector_size;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Cuda, TVectorSize > {
    constexpr static int value = 256 / TVectorSize;
  };
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int get_suggested_team_size( const int vector_size, Kokkos::OpenMP )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::OpenMP, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Serial )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Serial, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int get_suggested_team_size( const int vector_size, Kokkos::Threads )
  {
    return 1;
  }

  template< int TVectorSize >
  struct SuggestedTeamSize< Kokkos::Threads, TVectorSize > {
    constexpr static int value = 1;
  };
  #endif

  template< typename TExecSpace, int TTeamSize >
  struct SuggestedTeamWorkSize;

  template< typename TExecSpace >
  inline int get_team_work_size( const int team_size, TExecSpace )
  {
    throw std::runtime_error( "`get_team_work_size` is not implemented for requested execution space" );
  }

  #if defined(KOKKOS_ENABLE_CUDA)
  inline int get_team_work_size( const int team_size, Kokkos::Cuda )
  {
    return team_size;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Cuda, TTeamSize >
  {
    constexpr static int value = TTeamSize;
  };
  #endif

  #if defined(KOKKOS_ENABLE_OPENMP)
  inline int get_team_work_size( const int team_size, Kokkos::OpenMP )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::OpenMP, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif

  #if defined(KOKKOS_ENABLE_SERIAL)
  inline int get_team_work_size( const int team_size, Kokkos::Serial )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Serial, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif

  #if defined(KOKKOS_ENABLE_THREADS)
  inline int get_team_work_size( const int team_size, Kokkos::Threads )
  {
    return 16;
  }

  template< int TTeamSize >
  struct SuggestedTeamWorkSize< Kokkos::Threads, TTeamSize >
  {
    constexpr static int value = 16;
  };
  #endif
}  // namespace diverg
#endif  // DIVERG_RANGE_SPARSE_BASE_HPP_
