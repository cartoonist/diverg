/**
 *    @file  dindex_eti.hpp
 *   @brief  ETI (compiled-library) interface for distance index construction.
 *
 *  This is the Kokkos-free interface to the distance index construction that is
 *  compiled into `libdiverg`. It mirrors the header-only `dindex.hpp` (which
 *  includes the full templated implementation): downstream code that uses the
 *  ETI build includes *this* header instead of `dindex.hpp` and links against
 *  `libdiverg`, so it can build the adjacency, call `create_distance_index`,
 *  and assemble the result without pulling in Kokkos.
 *
 *  The `create_distance_index` is compiled with the default sparse
 *  configuration, so it runs on the device in a CUDA build. Its input and
 *  output are `CRSMatrix` host containers. Overloads are provided for the
 *  supported `<ordinal, size>` type combinations.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_DINDEX_ETI_HPP__
#define DIVERG_DINDEX_ETI_HPP__

#include <cstdint>

#include <diverg/crs_matrix.hpp>


namespace diverg {
  namespace util {
    /**
     *  @brief  Build the distance index of a region from its range adjacency.
     *
     *  @param  ra    The adjacency matrix of the region in range CRS format.
     *  @param  dmin  Lower bound of the valid distance range (inclusive).
     *  @param  dmax  Upper bound of the valid distance range (inclusive).
     *  @return The (uncompressed, mutable range) distance index.
     *
     *  NOTE: Uses the default sparse configuration; the execution space follows
     *  the Kokkos backend the library was built with.
     */
    CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t >
    create_distance_index(
        CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > const& ra,
        unsigned int dmin, unsigned int dmax );

    /** @overload for a 32-bit size type. */
    CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint32_t >
    create_distance_index(
        CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint32_t > const& ra,
        unsigned int dmin, unsigned int dmax );
  }  /* --- end of namespace util --- */
}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_DINDEX_ETI_HPP__ --- */
