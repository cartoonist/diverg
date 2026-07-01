/**
 *    @file  crs_matrix_device.hpp
 *   @brief  Kokkos view adapters for `CRSMatrix`.
 *
 *  CRSMatrix definition is separated into two header files: `crs_matrix.hpp`,
 *  which defines `CRSMatrix` as a Kokkos-free *host* container plus query and
 *  serialisation on Host, and this header file, which adds the *device* side.
 *
 *  This header file creates Kokkos views over a CRSMatrix's host storage and
 *  mirrors them to a device execution space. These are provided as functions
 *  so that downstream host translation units can include `crs_matrix.hpp`
 *  without pulling in Kokkos, while DiVerG internals (and the compiled ETI
 *  library) include this header to operate on matrices on device.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_CRS_MATRIX_DEVICE_HPP__
#define DIVERG_CRS_MATRIX_DEVICE_HPP__

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <Kokkos_Core.hpp>

#include "crs_matrix.hpp"


namespace diverg {
  namespace crs_matrix {
    /**
     *  @brief  Whether a host view can be created over a spec's `entries` storage.
     *
     *  Only in-memory contiguous storage is viewable; out-of-core (buffered) and
     *  compressed specialisations are not.
     */
    template< typename TSpec > struct is_entries_viewable : std::false_type { };
    template< > struct is_entries_viewable< Dynamic > : std::true_type { };
    template< > struct is_entries_viewable< RangeDynamic > : std::true_type { };

    /**
     *  @brief  Whether a host view can be created over a spec's `rowmap` storage.
     */
    template< typename TSpec > struct is_rowmap_viewable : std::false_type { };
    template< > struct is_rowmap_viewable< Dynamic > : std::true_type { };
    template< > struct is_rowmap_viewable< Buffered > : std::true_type { };
    template< > struct is_rowmap_viewable< RangeDynamic > : std::true_type { };
    template< > struct is_rowmap_viewable< RangeBuffered > : std::true_type { };

    /**
     *  @brief  Count non-zero values of a Range CRS matrix from device views.
     *
     *  Device counterpart of the host `nnz( entries, rowmap, RangeGroup )` in
     *  `crs_matrix.hpp`: sums each `[begin, end]` range length in parallel.
     *
     *  NOTE: Kept for reference. It is currently unused as the range sparse
     *  algorithms track and return the result nnz themselves.
     */
    template< typename TRowMapDeviceView, typename ...TArgs >
    inline typename TRowMapDeviceView::non_const_value_type
    nnz( Kokkos::View< TArgs... > d_entries, TRowMapDeviceView, RangeGroup /*tag*/ )
    {
      typedef Kokkos::View< TArgs... > entries_t;
      typedef TRowMapDeviceView row_map_t;
      typedef typename row_map_t::non_const_value_type size_type;
      typedef typename entries_t::execution_space execution_space;
      typedef Kokkos::RangePolicy< execution_space > policy_type;

      /* just check if both `Kokkos::View`s are on the same space */
      static_assert( std::is_same< typename entries_t::memory_space,
                                   typename row_map_t::memory_space >::value,
                     "different memory space" );

      auto ent_size = d_entries.extent( 0 );
      assert( ent_size % 2 == 0 );

      size_type nnz_counter = 0;
      Kokkos::parallel_reduce(
          "diverg::crs_matrix::nnz_range", policy_type( 0, ent_size / 2 ),
          KOKKOS_LAMBDA ( const uint64_t ii, size_type& nnz_local ) {
            nnz_local += d_entries( ii*2 + 1 ) - d_entries( ii*2 ) + 1;
          },
          nnz_counter );  // reduce to scalar is blocking

      return nnz_counter;
    }
  }  /* --- end of namespace crs_matrix --- */

  /* === HOST-SPACE UNMANAGED VIEW TYPES (over the matrix storage) === */
  template< typename TMatrix, typename TDeviceSpace >
  using entries_view_type = Kokkos::View<
      typename TMatrix::ordinal_type*, typename TDeviceSpace::array_layout,
      Kokkos::DefaultHostExecutionSpace, Kokkos::MemoryUnmanaged >;

  template< typename TMatrix, typename TDeviceSpace >
  using const_entries_view_type = Kokkos::View<
      const typename TMatrix::ordinal_type*, typename TDeviceSpace::array_layout,
      Kokkos::DefaultHostExecutionSpace, Kokkos::MemoryUnmanaged >;

  template< typename TMatrix, typename TDeviceSpace >
  using rowmap_view_type = Kokkos::View<
      typename TMatrix::size_type*, typename TDeviceSpace::array_layout,
      Kokkos::DefaultHostExecutionSpace, Kokkos::MemoryUnmanaged >;

  template< typename TMatrix, typename TDeviceSpace >
  using const_rowmap_view_type = Kokkos::View<
      const typename TMatrix::size_type*, typename TDeviceSpace::array_layout,
      Kokkos::DefaultHostExecutionSpace, Kokkos::MemoryUnmanaged >;

  /* === DEVICE-SPACE VIEW TYPES (managed mirrors) === */
  template< typename TMatrix, typename TDeviceSpace >
  using entries_device_view_type
      = Kokkos::View< typename TMatrix::ordinal_type*, TDeviceSpace >;

  template< typename TMatrix, typename TDeviceSpace >
  using const_entries_device_view_type
      = Kokkos::View< const typename TMatrix::ordinal_type*, TDeviceSpace >;

  template< typename TMatrix, typename TDeviceSpace >
  using rowmap_device_view_type
      = Kokkos::View< typename TMatrix::size_type*, TDeviceSpace >;

  template< typename TMatrix, typename TDeviceSpace >
  using const_rowmap_device_view_type
      = Kokkos::View< const typename TMatrix::size_type*, TDeviceSpace >;

  /* === HOST VIEWS OVER STORAGE === */
  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline auto
  entries_view( TMatrix& matrix, TDeviceSpace = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    static_assert( crs_matrix::is_entries_viewable< typename matrix_type::spec_type >::value,
                   "cannot create entries view for this CRSMatrix specialisation" );
    auto& e = matrix.get_entries();
    if ( e.size() != 0 )
      return entries_view_type< matrix_type, TDeviceSpace >( e.data(), e.size() );
    else
      return entries_view_type< matrix_type, TDeviceSpace >( nullptr, 0 );
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline auto
  entries_view( TMatrix const& matrix, TDeviceSpace = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    static_assert( crs_matrix::is_entries_viewable< typename matrix_type::spec_type >::value,
                   "cannot create entries view for this CRSMatrix specialisation" );
    auto const& e = matrix.get_entries();
    if ( e.size() != 0 )
      return const_entries_view_type< matrix_type, TDeviceSpace >( e.data(), e.size() );
    else
      return const_entries_view_type< matrix_type, TDeviceSpace >( nullptr, 0 );
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline auto
  rowmap_view( TMatrix& matrix, TDeviceSpace = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    static_assert( crs_matrix::is_rowmap_viewable< typename matrix_type::spec_type >::value,
                   "cannot create rowmap view for this CRSMatrix specialisation" );
    auto& r = matrix.get_rowmap();
    if ( r.size() != 0 )
      return rowmap_view_type< matrix_type, TDeviceSpace >( r.data(), r.size() );
    else
      return rowmap_view_type< matrix_type, TDeviceSpace >( nullptr, 0 );
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline auto
  rowmap_view( TMatrix const& matrix, TDeviceSpace = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    static_assert( crs_matrix::is_rowmap_viewable< typename matrix_type::spec_type >::value,
                   "cannot create rowmap view for this CRSMatrix specialisation" );
    auto const& r = matrix.get_rowmap();
    if ( r.size() != 0 )
      return const_rowmap_view_type< matrix_type, TDeviceSpace >( r.data(), r.size() );
    else
      return const_rowmap_view_type< matrix_type, TDeviceSpace >( nullptr, 0 );
  }

  /* === DEVICE MIRRORS === */
  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline entries_device_view_type< std::decay_t< TMatrix >, TDeviceSpace >
  entries_device_view( TMatrix& matrix, TDeviceSpace space = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    if ( matrix.get_entries().size() != 0 )
      return Kokkos::create_mirror_view_and_copy(
          TDeviceSpace{}, entries_view( matrix, space ) );
    else
      return entries_device_view_type< matrix_type, TDeviceSpace >{};
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline const_entries_device_view_type< std::decay_t< TMatrix >, TDeviceSpace >
  entries_device_view( TMatrix const& matrix, TDeviceSpace space = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    if ( matrix.get_entries().size() != 0 )
      return Kokkos::create_mirror_view_and_copy(
          TDeviceSpace{}, entries_view( matrix, space ) );
    else
      return const_entries_device_view_type< matrix_type, TDeviceSpace >{};
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline rowmap_device_view_type< std::decay_t< TMatrix >, TDeviceSpace >
  rowmap_device_view( TMatrix& matrix, TDeviceSpace space = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    if ( matrix.get_rowmap().size() != 0 )
      return Kokkos::create_mirror_view_and_copy(
          TDeviceSpace{}, rowmap_view( matrix, space ) );
    else
      return rowmap_device_view_type< matrix_type, TDeviceSpace >{};
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline const_rowmap_device_view_type< std::decay_t< TMatrix >, TDeviceSpace >
  rowmap_device_view( TMatrix const& matrix, TDeviceSpace space = {} )
  {
    typedef std::decay_t< TMatrix > matrix_type;
    if ( matrix.get_rowmap().size() != 0 )
      return Kokkos::create_mirror_view_and_copy(
          TDeviceSpace{}, rowmap_view( matrix, space ) );
    else
      return const_rowmap_device_view_type< matrix_type, TDeviceSpace >{};
  }

  /* === EMPTY DEVICE VIEWS (factory by matrix type) === */
  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline entries_device_view_type< TMatrix, TDeviceSpace >
  make_entries_device_view( TDeviceSpace = {} )
  {
    return entries_device_view_type< TMatrix, TDeviceSpace >{};
  }

  template< typename TMatrix, typename TDeviceSpace = Kokkos::DefaultExecutionSpace >
  inline rowmap_device_view_type< TMatrix, TDeviceSpace >
  make_rowmap_device_view( TDeviceSpace = {} )
  {
    return rowmap_device_view_type< TMatrix, TDeviceSpace >{};
  }

  /**
   *  @brief  Construct a (host) CRSMatrix from device views.
   *
   *  Mirrors the device `entries`/`rowmap` views into freshly-allocated host
   *  storage and constructs the matrix from it. Replaces the device-view
   *  constructors that previously lived on `CRSMatrix` (which required Kokkos).
   *
   *  NOTE: `TMatrix` must have viewable (in-memory contiguous) storage, i.e.
   *  `Dynamic`/`RangeDynamic`.
   */
  template< typename TMatrix, typename TEntriesDView, typename TRowmapDView >
  inline TMatrix
  make_crs_matrix( typename TMatrix::ordinal_type ncols, TEntriesDView d_entries,
                   TRowmapDView d_rowmap, typename TMatrix::size_type nnz )
  {
    typedef Kokkos::DefaultExecutionSpace device_space;
    static_assert( crs_matrix::is_entries_viewable< typename TMatrix::spec_type >::value
                   && crs_matrix::is_rowmap_viewable< typename TMatrix::spec_type >::value,
                   "make_crs_matrix requires a CRSMatrix with viewable storage" );

    typename TMatrix::entries_type h_entries = TMatrix::make_entries( d_entries.extent( 0 ) );
    typename TMatrix::rowmap_type h_rowmap = TMatrix::make_rowmap( d_rowmap.extent( 0 ) );

    entries_view_type< TMatrix, device_space > h_entries_view(
        h_entries.data(), h_entries.size() );
    rowmap_view_type< TMatrix, device_space > h_rowmap_view(
        h_rowmap.data(), h_rowmap.size() );
    // FIXME: An extra copy takes place here when the input views are already on
    // host (i.e. a host execution space): the deep copy is then host-to-host.
    // Probably, a specialised CRSMatrix with an internal `Kokkos::StaticCrsGraph`
    // will resolve the issue.
    Kokkos::deep_copy( h_entries_view, d_entries );
    Kokkos::deep_copy( h_rowmap_view, d_rowmap );

    return TMatrix( ncols, std::move( h_entries ), std::move( h_rowmap ), nnz );
  }
}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_CRS_MATRIX_DEVICE_HPP__ --- */
