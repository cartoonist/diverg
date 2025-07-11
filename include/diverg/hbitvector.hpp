/**
 *    @file  hbitvector.hpp
 *   @brief  HBitVector class definition header file.
 *
 *  This header file contains `HBitVector` class definition and related
 *  interface functions and type definitions.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Tue Jun 20, 2023  22:23
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2023, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_HBITVECTOR_HPP_
#define DIVERG_HBITVECTOR_HPP_

#include <climits>
#include <cstdint>
#include <limits>

#include <gum/basic_types.hpp>
#include <Kokkos_Core.hpp>

#include "utils.hpp"
#include "bits.hpp"
#include "range_sparse_base.hpp"


namespace diverg {
  template< typename TDevice >
  struct HBitVectorTraits {
    using size_type = uint32_t;
    using bitset_type = uint64_t;
  };

#if defined( KOKKOS_ENABLE_CUDA )
  template< >
  struct HBitVectorTraits< Kokkos::Cuda > {
    using size_type = uint32_t;
    using bitset_type = uint32_t;
  };
#endif

  struct Safe {};
  struct UnsafeAtBoundaries {};
  struct Unsafe {};

  template< typename TSpec >
  struct HBVThreadAccess {};

  template< typename TPartition >
  struct GetHBVThreadAccess;

  template< typename TSpec >
  struct GetHBVThreadAccess< ExecPartition< TSpec > > {
    using type = HBVThreadAccess< Unsafe >;
  };

  template<>
  struct GetHBVThreadAccess< ThreadSequentialPartition > {
    using type = HBVThreadAccess< Safe >;
  };

  template<>
  struct GetHBVThreadAccess< TeamSequentialPartition > {
    using type = HBVThreadAccess< UnsafeAtBoundaries >;
  };

  template< typename TPartition >
  using GetHBVThreadAccessType = typename GetHBVThreadAccess< TPartition >::type;

  struct TeamLevel {};
  struct ThreadLevel {};
  struct VectorLevel {};

  template< typename TSpec >
  struct HBVAccessLevel {};

  template< typename TPartition >
  struct GetHBVAccessLevel;

  template< typename TSpec >
  struct GetHBVAccessLevel< ExecPartition< TSpec > > {
    using type = HBVAccessLevel< TeamLevel >;
  };

  template<>
  struct GetHBVAccessLevel< ThreadSequentialPartition > {
    using type = HBVAccessLevel< ThreadLevel >;
  };

  template< typename TPartition >
  using GetHBVAccessLevelType = typename GetHBVAccessLevel< TPartition >::type;

  /**
   *  @brief  Banded Hierarchical (two-level) bit vector
   *
   *  @tparam  TL1_Size Size of the first level bit vector (faster)
   *
   *  NOTE: This class implements *banded* hbitvector centring around the given
   *  'centre' index. The specified size when constructing is the size of the
   *  band that includes L1 rather than the whole bitvector. The hbitvector
   *  only checks whether the accessed bits are within the band and there is no
   *  check or knowledge for the actual size of the vector.
   *
   *  The arrangement of two bit arrays for each level is shown below:
   *
   *  n = vector size, N = aligned vector size (multiple of bitset width)
   *  g = global index, b = L1 begin position, r = relative index, l = local index
   *  g: b b+1  ...  b+|L1|-1   b+|L1| ... n 0 1 ...  b-1
   *     | |               |    |          | | |       |
   *   [     L1 region      ] [        L2 region         ]
   *     | |               |    | | |                  |
   *  l: 0 1    ...    |L1|-1   0 1 2      ...     |L2|-1
   *
   *  r = ( ( N + g - b ) % N )
   *  l = r < |L1| ? r : r - |L1|
   */
  template< unsigned int TL1_Size = 2048, /*bits*/
            typename TDevice = Kokkos::DefaultExecutionSpace,
            typename TTrait = HBitVectorTraits< TDevice > >
  class HBitVector {
    public:
      /* === MEMBER TYPES === */
      using device_type = TDevice;
      using trait_type = TTrait;
      using execution_space = typename device_type::execution_space;
      using scratch_space = typename execution_space::scratch_memory_space;
      using policy_type = typename Kokkos::TeamPolicy< execution_space >;
      using member_type = typename policy_type::member_type;
      using size_type = typename trait_type::size_type;
      using bitset_type = typename trait_type::bitset_type;
      using view_type = bitset_type*;
      using value_type = bool;
      /* === MEMBER CONSTANTS === */
      static constexpr const size_type BITSET_WIDTH            = gum::widthof< bitset_type >::value;  // 64 (if uint64_t)
      static constexpr const unsigned short int BINDEX_SHIFT   = diverg::bits::hi( BITSET_WIDTH );    // 6  (if uint64_t)
      static constexpr const size_type BOFFSET_MASK            = BITSET_WIDTH - 1u; // 0x0000001f (if size_type=uint32_t)
      static constexpr const size_type INDEX_ALIGN_MASK        = ~BOFFSET_MASK;     // 0xffffffe0 (if size_type=uint32_t)
      static constexpr const bitset_type BITSET_ALL_NIL        = 0;                   // 0x0000000000000000 (if uint64_t)
      static constexpr const bitset_type BITSET_ALL_SET        = ~BITSET_ALL_NIL;     // 0xffffffffffffffff (if uint64_t)
      static constexpr const bitset_type BITSET_ONE            = 1u;
      static constexpr const size_type L1_SIZE                 = TL1_Size;
      static constexpr const size_type L1_SIZE_BYTES           = TL1_Size / CHAR_BIT;
      static constexpr const size_type L1_NUM_BITSETS          = TL1_Size >> BINDEX_SHIFT;
      static constexpr const std::size_t value_alignment       = Kokkos::max( sizeof( bitset_type ), alignof( bitset_type ) );
      static constexpr const std::size_t space_alignment       = Kokkos::max( value_alignment, static_cast< size_t >( scratch_space::ALIGN ) );
      /* === STATIC ASSERTS === */
      // Accepting 64-bit L1 size requires spending extra time checking for corner cases
      static_assert( ( TL1_Size >= ( BITSET_WIDTH << 1 ) ), "L1 size should be at least twice larger than bitset width" );
      static_assert( ( diverg::bits::cnt( TL1_Size ) == 1 ), "L1 size should be a power of 2" );
      static_assert( ( diverg::bits::cnt( BITSET_WIDTH ) == 1 ), "Bitset width should be a power of 2" );
      static_assert( ( TL1_Size <= std::numeric_limits< size_type >::max() ), "L1 size cannot fit in size type" );
      /* === DATA MEMBERS === */
      //size_type m_size;         //!< Size of the bit vector
      size_type m_x_size;       //!< Allocated *BANDED* size in bits (always a multiple of `BITSET_WIDTH`)
      size_type m_num_bitsets;  //!< Total number of bitsets in *BANDED* form
      size_type l1_begin_bidx;  //!< Bitset index of the first bitset in L1 (inclusive)
      size_type l1_begin;       //!< Index of the first bit reside in L1 (inclusive)
      view_type l1_data;        //!< First level bit vector view
      view_type l2_data;        //!< Second level bit vector view
    private:
      /**
       *  NOTE: For banded hbitvector, `n` should be band size.
       */
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n )
          : /*m_size( n ),*/ l1_begin_bidx( 0 ), l1_begin( 0 ),
            l1_data(), l2_data()
      {
        this->m_x_size = HBitVector::aligned_size( n );
        this->m_num_bitsets = HBitVector::bindex( this->m_x_size );
      }

      KOKKOS_INLINE_FUNCTION void
      allocate_data( const member_type& tm, HBVAccessLevel< TeamLevel > )
      {
        auto l1size = HBitVector::l1_scratch_size();
        this->l1_data = ( view_type )
          ( tm.team_scratch( 0 ).get_shmem_aligned( l1size, space_alignment ) );
        auto l2size = this->l2_scratch_size();
        if ( l2size != 0 ) {
          this->l2_data = ( view_type )
            ( tm.team_scratch( 1 ).get_shmem_aligned( l2size, space_alignment ) );
        }
      }

      KOKKOS_INLINE_FUNCTION void
      allocate_data( const member_type& tm, HBVAccessLevel< ThreadLevel > )
      {
        auto l1size = HBitVector::l1_scratch_size();
        this->l1_data = ( view_type )
          ( tm.thread_scratch( 0 ).get_shmem_aligned( l1size, space_alignment ) );
        auto l2size = this->l2_scratch_size();
        if ( l2size != 0 ) {
          this->l2_data = ( view_type )
            ( tm.thread_scratch( 1 ).get_shmem_aligned( l2size, space_alignment ) );
        }
      }
    public:
      /* === LIFECYCLE === */
      /**
       *  @brief Constructor for NON-banded hbitvector
       *
       *  @param tm     Kokkos team member
       *  @param n      Actual size
       *  @param centre Center index
       *  @param tag    Access level tag
       */
      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n, size_type centre,
                  HBVAccessLevel< TSpec > access_level_tag )
          : HBitVector( tm, n )
      {
        this->allocate_data( tm, access_level_tag );
        this->set_l1_centre_at( centre, n );
      }

      /**
       *  @brief Constructor for banded hbitvector
       *
       *  @param tm     Kokkos team member
       *  @param dim    HBitvector dimensions: (banded size, actual size)
       *  @param centre Center index
       *  @param tag    Access level tag
       */
      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, std::pair< size_type, size_type > dim,
                  size_type centre, HBVAccessLevel< TSpec > access_level_tag )
          : HBitVector( tm, dim.first /* banded size */ )
      {
        this->allocate_data( tm, access_level_tag );
        this->set_l1_centre_at( centre, dim.second /* actual size */ );
      }

      /**
       *  NOTE: For banded hbitvector, `n` should be band size.
       */
      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n,
                  HBVAccessLevel< TSpec > access_level_tag )
          : HBitVector( tm, n )
      {
        this->allocate_data( tm, access_level_tag );
      }

      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n, size_type centre,
                  ExecPartition< TSpec > )
          : HBitVector( tm, n, centre,
                        GetHBVAccessLevelType< ExecPartition< TSpec > >{} )
      { }

      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, std::pair< size_type, size_type > dim,
                  size_type centre, ExecPartition< TSpec > )
          : HBitVector( tm, dim, centre,
                        GetHBVAccessLevelType< ExecPartition< TSpec > >{} )
      { }

      template< typename TSpec >
      KOKKOS_FUNCTION
      HBitVector( const member_type& tm, size_type n, ExecPartition< TSpec > )
          : HBitVector( tm, n,
                        GetHBVAccessLevelType< ExecPartition< TSpec > >{} )
      { }
      /* === STATIC MEMBERS === */
      static KOKKOS_INLINE_FUNCTION size_type
      compute_l1_begin_bidx( size_type centre, size_type actual_size )
      {
        auto x_size = HBitVector::aligned_size( actual_size );
        auto num_bitsets = HBitVector::bindex( x_size );

        assert( centre < x_size );

        auto ctr_bidx = HBitVector::bindex( centre );
        // The begin index is set such that the L1 range would be (inclusive):
        //   [ctr_bidx-(L1_NUM_BITSETS/2)+1...ctr_bidx+(L1_NUM_BITSETS/2)]
        if ( L1_NUM_BITSETS < num_bitsets ) {
          // lflank: L1 left flank size relative to centre
          auto lflank = ( L1_NUM_BITSETS >> 1 ) - 1;
          // rfit_bidx: right-most bitset index fitting the whole L1
          auto rfit_bidx = num_bitsets - L1_NUM_BITSETS;
          auto pb_bidx = ( ( ctr_bidx > lflank ) ? ( ctr_bidx - lflank ) : 0 );
          // for values of `centre` being closer to the end, l1 covers the last `TL1_Size` bits
          return Kokkos::min( rfit_bidx, pb_bidx );
        }
        return 0;
      }

      /**
       *   @brief Get bitset index of the bit at `idx`
       *
       *   @param idx Bit index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      bindex( size_type idx ) noexcept
      {
        return idx >> BINDEX_SHIFT;
      }

      /**
       *   @brief Get bitset offset of the bit at `idx`
       *
       *   @param idx Bit index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      boffset( size_type idx ) noexcept
      {
        return idx & BOFFSET_MASK;
      }

      /**
       *   @brief Get index of the starting bit of the bitset with index `bidx`
       *
       *   @param bidx Bitset index
       */
      static KOKKOS_INLINE_FUNCTION size_type
      start_index( size_type bidx ) noexcept
      {
        return bidx << BINDEX_SHIFT;
      }

      /**
       *   @brief Get the left closest 'aligned index' relative to index `idx`
       *
       *   @param idx Bit index
       *
       *   Aligned index of an index, lets say `idx`, is the index of the
       *   starting bit of the bitset which includes the bit at `idx`.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_index( size_type idx ) noexcept
      {
        return idx & INDEX_ALIGN_MASK;
      }

      /**
       *   @brief Get the right closest 'aligned index' relative to index `idx`
       *
       *   @param idx Bit index
       *
       *   @return the index of the starting bit of the bitset immediately
       *   after the one including the bit at `idx` unless `idx` is itself a
       *   starting bit of a bitset.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_index_ceil( size_type idx ) noexcept
      {
        return ( idx + ( BITSET_WIDTH - 1 ) ) & ( INDEX_ALIGN_MASK );
      }

      static KOKKOS_INLINE_FUNCTION size_type
      space_aligned_size( size_type n ) noexcept
      {
        // 0x00000007 if `space_alignment==8` and `size_type=uint32_t`
        size_type hi_mask = HBitVector::space_alignment - 1;
        // 0xfffffff8 if `space_alignment==8` and `size_type=uint32_t`
        size_type lo_mask = ~hi_mask;
        // the result would be `n` if it is already a multiple of
        // `space_alignment`; otherwise the next closest one.
        return ( n + hi_mask ) & lo_mask;
      }

      /**
       *   @brief Return aligned size for a vector of `n` bits
       *
       *   The aligned size is the smallest multiple of `BITSET_WIDTH` which is
       *   larger than the actual size of the bitvector. The minimum aligned
       *   size is the size of L1.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      aligned_size( size_type n )
      {
        auto asize = ( n + ( BITSET_WIDTH - 1 ) ) & ( INDEX_ALIGN_MASK );  // aligned_index_ceil( n );
        return Kokkos::max( asize, TL1_Size );
      }

      /**
       *   @brief Return the number of allocated bitsets for L1
       *
       *   NOTE: The vector itself might occupy less bitsets than the number of
       *   allocated ones as the L1 allocated size is fixed in compile time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_num_bitsets( ) noexcept
      {
        return L1_NUM_BITSETS;
      }

      /**
       *   @brief Return the allocated size of L1 in bits
       *
       *   NOTE: The vector itself might occupy less bits than the value return
       *   by this function as the L1 allocated size is fixed in compile time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_size( ) noexcept
      {
        return TL1_Size;
      }

      /**
       *   @brief Return the allocated size of L1 in bytes
       *
       *   NOTE: The vector itself might occupy less bytes than the value
       *   return by this function as the L1 allocated size is fixed in compile
       *   time.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l1_scratch_size( ) noexcept
      {
        return L1_SIZE_BYTES;
      }

      /**
       *   @brief Return the required number of bitsets for a vector of `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      num_bitsets( size_type n ) noexcept
      {
        auto x_size = HBitVector::aligned_size( n );
        return HBitVector::bindex( x_size );
      }

      /**
       *   @brief Return the required number of bitsets in L2 for a vector of
       *          `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_num_bitsets( size_type n ) noexcept
      {
        auto nbitsets = HBitVector::num_bitsets( n );
        return ( nbitsets > HBitVector::l1_num_bitsets() )
                   ? nbitsets - HBitVector::l1_num_bitsets()
                   : 0;
      }

      /**
       *   @brief Return the required L2 size in bits for a vector of `n` bits
       *
       *   NOTE: The vector itself might occupy less bits in L2 than the value
       *   return by this function which is the size need to be *allocated*.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_size( size_type n ) noexcept
      {
        auto x_size = HBitVector::aligned_size( n );
        return ( x_size > HBitVector::l1_size() )
                   ? x_size - HBitVector::l1_size()
                   : 0;
      }

      /**
       *   @brief Return the required L2 size in bytes for a vector of `n` bits
       *
       *   NOTE: The vector itself might occupy less bytes in L2 than the value
       *   return by this function which is the size need to be *allocated*.
       */
      static KOKKOS_INLINE_FUNCTION size_type
      l2_scratch_size( size_type n ) noexcept
      {
        return HBitVector::l2_size( n ) / CHAR_BIT;
      }

      /**
       *   @brief Return the required number of bytes for a vector of `n` bits
       */
      static KOKKOS_INLINE_FUNCTION size_type
      capacity( size_type n ) noexcept
      {
        return HBitVector::aligned_size( n ) / CHAR_BIT;
      }

      template< typename TPolicy >
      static inline TPolicy
      set_scratch_size( TPolicy& policy, size_type n,
                        HBVAccessLevel< TeamLevel > )
      {
        auto l1size = HBitVector::l1_scratch_size();
        auto l2size_aln = HBitVector::space_aligned_size(
            HBitVector::l2_scratch_size( n ) );
        policy.set_scratch_size( 0, Kokkos::PerTeam( l1size ) );
        if ( l2size_aln != 0 ) {
          policy.set_scratch_size( 1, Kokkos::PerTeam( l2size_aln ) );
        }
        return policy;
      }

      template< typename TPolicy >
      static inline TPolicy
      set_scratch_size( TPolicy& policy, size_type n,
                        HBVAccessLevel< ThreadLevel > )
      {
        auto l1size = HBitVector::l1_scratch_size();
        auto l2size_aln = HBitVector::space_aligned_size(
            HBitVector::l2_scratch_size( n ) );
        policy.set_scratch_size( 0, Kokkos::PerThread( l1size ) );
        if ( l2size_aln != 0 ) {
          policy.set_scratch_size( 1, Kokkos::PerThread( l2size_aln ) );
        }
        return policy;
      }

      template< typename TPolicy, typename TSpec >
      static inline TPolicy
      set_scratch_size( TPolicy& policy, size_type n, ExecPartition< TSpec > )
      {
        return HBitVector::set_scratch_size(
            policy, n, GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      static KOKKOS_INLINE_FUNCTION size_type
      cnt( bitset_type x ) noexcept
      {
        return Kokkos::Experimental::popcount_builtin( x );  // Requires Kokkos >=4.1.00
      }

      /**
       *   @brief Calculate the position of the i-th rightmost 1 bit in `x`
       *
       *   @param x Input bitset value
       *   @param i Argument i must be in the range [1..cnt(x)].
       */
      static KOKKOS_INLINE_FUNCTION size_type
      sel( bitset_type x, size_type i )
      {
        KOKKOS_IF_ON_DEVICE( ( return HBitVector::sel_device( x, i ); ) )
        KOKKOS_IF_ON_HOST( ( return diverg::bits::sel( x, i ); ) )
      }

#if defined( KOKKOS_ENABLE_CUDA )
      /**
       *   @brief Calculate the position of the i-th rightmost 1 bit in `x`
       *          (on CUDA)
       *
       *   @param x Input bitset value
       *   @param i Argument i must be in the range [1..cnt(x)].
       */
      static KOKKOS_INLINE_FUNCTION size_type
      sel_device( bitset_type x, size_type i )
      {
        if constexpr ( BITSET_WIDTH == 64 ) {
          uint32_t lsw = x & 0xffffffff;
          auto cnt_lsw = HBitVector::cnt( lsw );
          if ( i <= cnt_lsw ) return __fns( lsw, 0, i );
          else {
            uint32_t msw = ( x >> 32 ) & 0xffffffff;
            return __fns( msw, 0, i - cnt_lsw );
          }
        }
        else {
          return __fns( x, 0, i );
        }
        DIVERG_BUILTIN_UNREACHABLE();
      }
#endif

      static KOKKOS_INLINE_FUNCTION bitset_type
      msb( bitset_type x ) noexcept
      {
        return x >> ( BITSET_WIDTH - 1 );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      lsb( bitset_type x ) noexcept
      {
        return x % 2;
      }

      // Following methods are taken from sdsl::bits
      static KOKKOS_INLINE_FUNCTION size_type
      cnt10( bitset_type x, bitset_type c ) noexcept
      {
        return HBitVector::cnt( ( ( x << 1 ) | c ) & ( ~x ) );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      map10( bitset_type x, bitset_type c ) noexcept
      {
        return ( ( ( x << 1 ) | c ) & ( ~x ) );
      }

      static KOKKOS_INLINE_FUNCTION size_type
      cnt01( bitset_type x, bitset_type c ) noexcept
      {
        return HBitVector::cnt( ( x ^ ( ( x << 1 ) | c ) ) & x );
      }

      static KOKKOS_INLINE_FUNCTION bitset_type
      map01( bitset_type x, bitset_type c ) noexcept
      {
        return ( ( x ^ ( ( x << 1 ) | c ) ) & x );
      }
      /* === OPERATORS === */
      /**
       *   @brief Get bitset by *global* bitset index.
       */
      KOKKOS_INLINE_FUNCTION bitset_type&
      operator()( size_type bidx ) const
      {
        auto r_bidx = this->relative_bitset( bidx );
        if ( r_bidx < HBitVector::l1_num_bitsets() ) {
          return this->l1_data[ r_bidx ];
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          return this->l2_data[ r_bidx ];
        }
      }

      KOKKOS_INLINE_FUNCTION value_type
      operator[]( size_type idx ) const
      {
        auto r_idx = this->relative_idx( idx );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto offset = HBitVector::boffset( r_idx );
        if ( r_idx < HBitVector::l1_size() ) {
          return ( this->l1_data[ r_bidx ] >> offset ) /*& BITSET_ONE*/;
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          return ( this->l2_data[ r_bidx ] >> offset ) /*& BITSET_ONE*/;
        }
      }
      /* === ACCESSORS === */
      /*
      KOKKOS_INLINE_FUNCTION size_type
      size( ) const
      {
        return this->m_size;
      }
      */

      /**
       *   @brief Return aligned size of the bitvector
       *
       *   The aligned size is the smallest multiple of `BITSET_WIDTH` which is
       *   larger than the actual size of the bitvector.
       */
      KOKKOS_INLINE_FUNCTION size_type
      aligned_size( ) const
      {
        return this->m_x_size;
      }

      KOKKOS_INLINE_FUNCTION size_type
      num_bitsets( ) const
      {
        return this->m_num_bitsets;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l1_begin_bindex( ) const
      {
        return this->l1_begin_bidx;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l1_begin_idx( ) const
      {
        return this->l1_begin;
      }
      /* === MUTATORS === */
      KOKKOS_INLINE_FUNCTION void
      set_l1_centre_at( size_type centre, size_type actual_size )
      {
        this->l1_begin_bidx
            = HBitVector::compute_l1_begin_bidx( centre, actual_size );
        this->l1_begin = HBitVector::start_index( this->l1_begin_bidx );
      }

      KOKKOS_INLINE_FUNCTION void
      set_l1_at( size_type idx )
      {
        this->l1_begin = HBitVector::aligned_index( idx );
        this->l1_begin_bidx = HBitVector::bindex( idx );
      }
      /* === METHODS === */
      /**
       *   @brief Returns the allocated size (L1+L2) in bytes
       */
      KOKKOS_INLINE_FUNCTION size_type
      capacity( ) const noexcept
      {
        return this->m_x_size / CHAR_BIT;
      }

      KOKKOS_INLINE_FUNCTION size_type
      l2_num_bitsets( ) const noexcept
      {
        auto nbitsets = this->num_bitsets();
        return ( nbitsets > HBitVector::l1_num_bitsets() )
                   ? nbitsets - HBitVector::l1_num_bitsets()
                   : 0;
      }

      /**
       *   @brief Return the allocated size of L2 in bits
       *
       *   NOTE: The vector itself might occupy less bits in L2 than the value
       *   return by this function which is the *allocated* size.
       */
      KOKKOS_INLINE_FUNCTION size_type
      l2_size( ) const noexcept
      {
        return ( this->m_x_size > HBitVector::l1_size() )
                   ? this->m_x_size - HBitVector::l1_size()
                   : 0;
      }

      /**
       *   @brief Return the allocated size of L2 in bytes
       *
       *   NOTE: The vector itself might occupy less bytes in L2 than the value
       *   return by this function which is the *allocated* size.
       */
      KOKKOS_INLINE_FUNCTION size_type
      l2_scratch_size( ) const noexcept
      {
        return this->l2_size() / CHAR_BIT;
      }

      /**
       *   @brief Return the relative index of `idx`
       *
       *   NOTE: The bitvector leaves the unused bit space [size, aln_size)
       *   untouched meaning that rearranging the bits for multi-level
       *   decompistion is done as if all allocated bits are used.
       */
      KOKKOS_INLINE_FUNCTION size_type
      relative_idx( size_type idx ) const noexcept
      {
        // NOTE: Using if-else is faster than computing modulo on both CPU and GPU
        //return ( this->m_x_size + idx - this->l1_begin ) % this->m_x_size;
        if ( this->l1_begin <= idx ) return idx - this->l1_begin;
        else return this->m_x_size + idx - this->l1_begin;
      }

      KOKKOS_INLINE_FUNCTION size_type
      relative_bitset( size_type bidx ) const noexcept
      {
        // NOTE: Using if-else is faster than computing modulo on both CPU and GPU
        //return ( this->m_num_bitsets + bidx - this->l1_begin_bidx )
        //       % this->m_num_bitsets;
        if ( this->l1_begin_bidx <= bidx ) return bidx - this->l1_begin_bidx;
        else return this->m_num_bitsets + bidx - this->l1_begin_bidx;
      }

      /**
       *   @brief Zero all bitsets in L1 (team-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l1( const member_type& tm, HBVAccessLevel< TeamLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, HBitVector::l1_num_bitsets() ),
            [=]( const uint64_t j ) { this->l1_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero all bitsets in L1 (thread-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l1( const member_type& tm, HBVAccessLevel< ThreadLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange( tm, HBitVector::l1_num_bitsets() ),
            [=]( const uint64_t j ) { this->l1_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero all bitsets in L1 (auto-detect by partition)
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      clear_l1( const member_type& tm, ExecPartition< TSpec > ) noexcept
      {
        this->clear_l1( tm, GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Zero all bitsets in L2 (team-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, HBVAccessLevel< TeamLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, this->l2_num_bitsets() ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero all bitsets in L2 (thread-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, HBVAccessLevel< ThreadLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange( tm, this->l2_num_bitsets() ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero all bitsets in L2 (auto-detect by partition)
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, ExecPartition< TSpec > ) noexcept
      {
        this->clear_l2( tm, GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Zero range of bitsets [ls_bidx, lf_bidx) on L2 (team-level)
       *
       *   @param tm  Team Policy member
       *   @param ls_bidx Local bitset index on L2 (inclusive)
       *   @param lf_bidx Local bitset index on L2 (exclusive)
       *
       *   NOTE: Both input indices are local bitset indices.
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, size_type ls_bidx,
                size_type lf_bidx, HBVAccessLevel< TeamLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [ls_bidx, lf_bidx) on L2 (thread-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, size_type ls_bidx,
                size_type lf_bidx, HBVAccessLevel< ThreadLevel > ) noexcept
      {
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [ls_bidx, lf_bidx) on L2
       *          (auto-detect by partition)
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      clear_l2( const member_type& tm, size_type ls_bidx,
                size_type lf_bidx, ExecPartition< TSpec > ) noexcept
      {
        this->clear_l2( tm, ls_bidx, lf_bidx,
                        GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Zero range of bitsets [s_bidx, f_bidx) which occurs on L2
       *          (team-level)
       *
       *   @param tm  Team Policy member
       *   @param s_bidx Global bitset index (inclusive)
       *   @param f_bidx Global bitset index (exclusive)
       *
       *   NOTE: Both input indices are global bitset indices (i.e. not local).
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_bidx( const member_type& tm, size_type s_bidx,
                        size_type f_bidx, HBVAccessLevel< TeamLevel > ) noexcept
      {
        assert( f_bidx != 0 );

        auto rs_bidx = this->relative_bitset( s_bidx );
        auto rf_bidx = this->relative_bitset( f_bidx - 1 ) + 1;
        auto ls_bidx = rs_bidx - HBitVector::l1_num_bitsets();
        auto lf_bidx = rf_bidx - HBitVector::l1_num_bitsets();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [s_bidx, f_bidx) which occurs on L2
       *          (thread-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_bidx( const member_type& tm, size_type s_bidx,
                        size_type f_bidx, HBVAccessLevel< ThreadLevel > ) noexcept
      {
        assert( f_bidx != 0 );

        auto rs_bidx = this->relative_bitset( s_bidx );
        auto rf_bidx = this->relative_bitset( f_bidx - 1 ) + 1;
        auto ls_bidx = rs_bidx - HBitVector::l1_num_bitsets();
        auto lf_bidx = rf_bidx - HBitVector::l1_num_bitsets();
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets [s_bidx, f_bidx) which occurs on L2
       *          (auto-detect by partition)
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_bidx( const member_type& tm, size_type s_bidx,
                        size_type f_bidx, ExecPartition< TSpec > ) noexcept
      {
        this->clear_l2_by_bidx( tm, s_bidx, f_bidx,
                                GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Zero range of bitsets on L2 covering bit range [s_idx, f_idx)
       *          (team-level)
       *
       *   @param tm  Team Policy member
       *   @param s_idx Global bit index in the vector (inclusive)
       *   @param f_idx Global bit index in the vector (exclusive)
       *
       *   NOTE: The function clears whole bitsets covering the range including
       *   the bitset on which `f_idx` lies unless `f_idx` is an aligned index.
       *
       *   NOTE: `s_idx` and `f_idx` are global *bit* indices (i.e. not local).
       *
       *   NOTE: Caller should make sure that the range does not span over L1.
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_idx( const member_type& tm, size_type s_idx,
                       size_type f_idx, HBVAccessLevel< TeamLevel > ) noexcept
      {
        assert( f_idx != 0 );

        auto rs_idx = this->relative_idx( s_idx );
        // `f_idx` can be equal to L1 begin index which would be mapped to 0
        // instead of |L1|+|L2| when computing the relative index.
        // `rel( f_idx - 1 ) + 1` gives the right answer without branching. The
        // following line combines this trick with computing the aligned index
        // ceiling method.
        auto rf_idx = ( this->relative_idx( f_idx - 1 ) + BITSET_WIDTH )
                      & ( INDEX_ALIGN_MASK );
        auto ls_idx = rs_idx - HBitVector::l1_size();
        auto lf_idx = rf_idx - HBitVector::l1_size();
        auto ls_bidx = HBitVector::bindex( ls_idx );
        auto lf_bidx = HBitVector::bindex( lf_idx );
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets on L2 covering bit range [s_idx, f_idx)
       *          (thread-level)
       */
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_idx( const member_type& tm, size_type s_idx,
                       size_type f_idx, HBVAccessLevel< ThreadLevel > ) noexcept
      {
        assert( f_idx != 0 );

        auto rs_idx = this->relative_idx( s_idx );
        // `f_idx` can be equal to L1 begin index which would be mapped to 0
        // instead of |L1|+|L2| when computing the relative index.
        // `rel( f_idx - 1 ) + 1` gives the right answer without branching. The
        // following line combines this trick with computing the aligned index
        // ceiling method.
        auto rf_idx = ( this->relative_idx( f_idx - 1 ) + BITSET_WIDTH )
                      & ( INDEX_ALIGN_MASK );
        auto ls_idx = rs_idx - HBitVector::l1_size();
        auto lf_idx = rf_idx - HBitVector::l1_size();
        auto ls_bidx = HBitVector::bindex( ls_idx );
        auto lf_bidx = HBitVector::bindex( lf_idx );
        Kokkos::parallel_for(
            Kokkos::ThreadVectorRange( tm, ls_bidx, lf_bidx ),
            [=]( const uint64_t j ) { this->l2_data[ j ] = 0; } );
      }

      /**
       *   @brief Zero range of bitsets on L2 covering bit range [s_idx, f_idx)
       *          (auto-detect by partition)
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      clear_l2_by_idx( const member_type& tm, size_type s_idx,
                       size_type f_idx, ExecPartition< TSpec > ) noexcept
      {
        this->clear_l2_by_idx( tm, s_idx, f_idx,
                               GetHBVAccessLevelType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Set a bit in the vector with Safe access
       *          (e.g. in ThreadSequentialParititon)
       *
       *   NOTE: The function should be called by a single thread/lane once.
       *
       *   NOTE: All write access are NON-atomic.
       */
      KOKKOS_INLINE_FUNCTION void
      set( size_type idx, HBVThreadAccess< Safe > ) noexcept
      {
        auto r_idx = this->relative_idx( idx );
        assert( r_idx < this->m_x_size );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto mask = BITSET_ONE << HBitVector::boffset( r_idx );

        if ( r_idx < HBitVector::l1_size() ) {  // most probable
          this->l1_data[ r_bidx ] |= mask;
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          this->l2_data[ r_bidx ] |= mask;
        }
      }

      /**
       *   @brief Set a bit in the vector
       *          (TeamSequential/ThreadParallel/TeamFlatParallel/etc)
       *
       *   NOTE: The function should be called by a single thread to be run on
       *         a vector lane (sequential).
       *
       *   NOTE: All write access are atomic.
       */
      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION std::enable_if_t< !std::is_same< TSpec, Safe >::value >
      set( size_type idx, HBVThreadAccess< TSpec > ) noexcept
      {
        auto r_idx = this->relative_idx( idx );
        assert( r_idx < this->m_x_size );
        auto r_bidx = HBitVector::bindex( r_idx );
        auto mask = BITSET_ONE << HBitVector::boffset( r_idx );

        if ( r_idx < HBitVector::l1_size() ) {  // most probable
          Kokkos::atomic_or( &this->l1_data[ r_bidx ], mask );
        }
        else {
          r_bidx -= HBitVector::l1_num_bitsets();
          Kokkos::atomic_or( &this->l2_data[ r_bidx ], mask );
        }
      }

      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      set( size_type idx, ExecPartition< TSpec > ) noexcept
      {
        this->set( idx, GetHBVThreadAccessType< ExecPartition< TSpec > >{} );
      }

      /**
       *   @brief Set a range of bits in the vector with Safe access
       *          (e.g. in ThreadSequentialPartition)
       *
       *   NOTE: The function should be called by a team to be run by a single
       *         thread (vector parallelism).
       */
      KOKKOS_INLINE_FUNCTION void
      set( const member_type& tm, size_type s_idx, size_type f_idx,
            HBVThreadAccess< Safe > tag ) noexcept
      {
        assert( s_idx <= f_idx );

        if ( s_idx == f_idx ) {
          Kokkos::single( Kokkos::PerThread( tm ), [=]() {
            this->set( s_idx, tag );
          } );
          return;
        }

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx );

        assert( rf_idx < this->m_x_size );

        auto setbits =
          [&tm]( auto data_ptr, auto ls_idx, auto lf_idx ) {
            auto ls_bidx = HBitVector::bindex( ls_idx );
            auto lf_bidx = HBitVector::bindex( lf_idx );

            if ( ls_bidx != lf_bidx ) {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto mask = ( BITSET_ALL_SET << s_offset );
                data_ptr[ ls_bidx ] |= mask;
              } );

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange( tm, ls_bidx + 1, lf_bidx ),
                  [=]( const uint64_t k ) {
                    data_ptr[ k ] |= BITSET_ALL_SET;
                  } );

              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ALL_SET >> ( BITSET_WIDTH - 1u - f_offset ) );
                data_ptr[ lf_bidx ] |= mask;
              } );
            }
            else {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ONE << ( f_offset - s_offset ) );
                mask = ( ( ( mask << 1 ) - 1 ) << s_offset );
                data_ptr[ ls_bidx ] |= mask;
              } );
            }
          };

        if ( rs_idx < l1_size() && rf_idx < l1_size() ) {  // the range is in L1
          setbits( this->l1_data, rs_idx, rf_idx );
        }
        else if ( rs_idx < l1_size() && l1_size() <= rf_idx ) {
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, rs_idx, l1_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rs_idx <= rf_idx ) {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l2_data, ls_idx, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rf_idx < l1_size() ) {
          auto ls_idx = rs_idx - l1_size();
          setbits( this->l1_data, 0, rf_idx );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
        }
        else {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, 0, l1_size() - 1 );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
      }

      /**
       *   @brief Set a range of bits in the vector with UnsafeAtBoundaries access
       *          (e.g. in TeamSequentialPartition)
       *
       *   NOTE: The function should be called by a team to be run by a single
       *         thread (vector parallelism).
       *
       *   NOTE: This function assumes that no other thread setting any bits in
       *         range [s_idx, f_idx] except for bitsets at end points.
       */
      KOKKOS_INLINE_FUNCTION void
      set( const member_type& tm, size_type s_idx, size_type f_idx,
            HBVThreadAccess< UnsafeAtBoundaries > tag ) noexcept
      {
        assert( s_idx <= f_idx );

        if ( s_idx == f_idx ) {
          Kokkos::single( Kokkos::PerThread( tm ), [=]() {
            this->set( s_idx, tag );
          } );
          return;
        }

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx );

        assert( rf_idx < this->m_x_size );

        auto setbits =
          [&tm]( auto data_ptr, auto ls_idx, auto lf_idx ) {
            auto ls_bidx = HBitVector::bindex( ls_idx );
            auto lf_bidx = HBitVector::bindex( lf_idx );

            if ( ls_bidx != lf_bidx ) {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto mask = ( BITSET_ALL_SET << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange( tm, ls_bidx + 1, lf_bidx ),
                  [=]( const uint64_t k ) {
                    data_ptr[ k ] |= BITSET_ALL_SET;
                  } );

              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ALL_SET >> ( BITSET_WIDTH - 1u - f_offset ) );
                Kokkos::atomic_or( &data_ptr[ lf_bidx ], mask );
              } );
            }
            else {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ONE << ( f_offset - s_offset ) );
                mask = ( ( ( mask << 1 ) - 1 ) << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );
            }
          };

        if ( rs_idx < l1_size() && rf_idx < l1_size() ) {  // the range is in L1
          setbits( this->l1_data, rs_idx, rf_idx );
        }
        else if ( rs_idx < l1_size() && l1_size() <= rf_idx ) {
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, rs_idx, l1_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rs_idx <= rf_idx ) {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l2_data, ls_idx, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rf_idx < l1_size() ) {
          auto ls_idx = rs_idx - l1_size();
          setbits( this->l1_data, 0, rf_idx );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
        }
        else {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, 0, l1_size() - 1 );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
      }

      /**
       *   @brief Set a range of bits in the vector with Unsafe access
       *          (e.g. in ThreadParallel)
       *
       *   NOTE: The function should be called by a team to be run by a single
       *         thread (vector parallelism).
       */
      KOKKOS_INLINE_FUNCTION void
      set( const member_type& tm, size_type s_idx, size_type f_idx,
            HBVThreadAccess< Unsafe > tag ) noexcept
      {
        assert( s_idx <= f_idx );

        if ( s_idx == f_idx ) {
          Kokkos::single( Kokkos::PerThread( tm ), [=]() {
            this->set( s_idx, tag );
          } );
          return;
        }

        auto rs_idx = this->relative_idx( s_idx );
        auto rf_idx = this->relative_idx( f_idx );

        assert( rf_idx < this->m_x_size );

        auto setbits =
          [&tm]( auto data_ptr, auto ls_idx, auto lf_idx ) {
            auto ls_bidx = HBitVector::bindex( ls_idx );
            auto lf_bidx = HBitVector::bindex( lf_idx );

            if ( ls_bidx != lf_bidx ) {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto mask = ( BITSET_ALL_SET << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );

              Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange( tm, ls_bidx + 1, lf_bidx ),
                  [=]( const uint64_t k ) {
                    Kokkos::atomic_or( &data_ptr[ k ], BITSET_ALL_SET );
                  } );

              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ALL_SET >> ( BITSET_WIDTH - 1u - f_offset ) );
                Kokkos::atomic_or( &data_ptr[ lf_bidx ], mask );
              } );
            }
            else {
              Kokkos::single( Kokkos::PerThread( tm ), [=]() {
                auto s_offset = HBitVector::boffset( ls_idx );
                auto f_offset = HBitVector::boffset( lf_idx );
                auto mask = ( BITSET_ONE << ( f_offset - s_offset ) );
                mask = ( ( ( mask << 1 ) - 1 ) << s_offset );
                Kokkos::atomic_or( &data_ptr[ ls_bidx ], mask );
              } );
            }
          };

        if ( rs_idx < l1_size() && rf_idx < l1_size() ) {  // the range is in L1
          setbits( this->l1_data, rs_idx, rf_idx );
        }
        else if ( rs_idx < l1_size() && l1_size() <= rf_idx ) {
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, rs_idx, l1_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rs_idx <= rf_idx ) {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l2_data, ls_idx, lf_idx );
        }
        else if ( l1_size() <= rs_idx && rf_idx < l1_size() ) {
          auto ls_idx = rs_idx - l1_size();
          setbits( this->l1_data, 0, rf_idx );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
        }
        else {
          auto ls_idx = rs_idx - l1_size();
          auto lf_idx = rf_idx - l1_size();
          setbits( this->l1_data, 0, l1_size() - 1 );
          setbits( this->l2_data, ls_idx, this->l2_size() - 1 );
          setbits( this->l2_data, 0, lf_idx );
        }
      }

      template< typename TSpec >
      KOKKOS_INLINE_FUNCTION void
      set( const member_type& tm, size_type s_idx, size_type f_idx,
           ExecPartition< TSpec > ) noexcept
      {
        this->set( tm, s_idx, f_idx,
                    GetHBVThreadAccessType< ExecPartition< TSpec > >{} );
      }
  };
}  /* --- end of namespace diverg --- */


#endif // DIVERG_HBITVECTOR_HPP_
