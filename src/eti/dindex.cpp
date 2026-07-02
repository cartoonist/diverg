/**
 *    @file  src/eti/dindex.cpp
 *   @brief  Compiled definitions of the distance index construction (ETI).
 *
 *  This source file defines `diverg::util::create_distance_index` which is
 *  declared Kokkos-free in `dindex_eti.hpp` for the supported `<ordinal, size>`
 *  type combinations. It forwards the call to the templated implementation in
 *  `dindex.hpp` with the *default sparse configuration*.
 *  
 *  This translation unit pulls in Kokkos, so it is compiled by the device
 *  compiler when CUDA is enabled and lives in the compiled `libdiverg` library.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <diverg/dindex.hpp>
#include <diverg/dindex_eti.hpp>


namespace diverg {
  namespace util {
    CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t >
    create_distance_index(
        CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t > const& ra,
        unsigned int dmin, unsigned int dmax )
    {
      return create_distance_index( ra, dmin, dmax, DefaultSparseConfiguration{} );
    }

    CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint32_t >
    create_distance_index(
        CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint32_t > const& ra,
        unsigned int dmin, unsigned int dmax )
    {
      return create_distance_index( ra, dmin, dmax, DefaultSparseConfiguration{} );
    }
  }  /* --- end of namespace util --- */
}  /* --- end of namespace diverg --- */
