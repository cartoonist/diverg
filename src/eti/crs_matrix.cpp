/**
 *    @file  src/eti/crs_matrix.cpp
 *   @brief  Explicit template instantiations for CRSMatrix.
 *
 *  This file contains ETI definitions for CRSMatrix. CRSMatrix has no device
 *  functions, so it is independent of CUDA (can be compiled with g++/clang++).
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Tue Mar 31, 2026  00:00
 *  Organization:  Universität Bielefeld
 *     Copyright:  Copyright (c) 2026, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

/* NOTE: Define this guard BEFORE including `crs_matrix.hpp` so that it skips
    the extern declarations (this source file provides the definitions). */
#define DIVERG_CRSMATRIX_ETI_INST__

#include <diverg/crs_matrix.hpp>

namespace diverg {

  /* --- Basic group --- */
  template class CRSMatrix< crs_matrix::Dynamic,       bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::Buffered,      bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::FullyBuffered, bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::Compressed,    bool, uint32_t, uint64_t >;

  /* --- Range group --- */
  template class CRSMatrix< crs_matrix::RangeDynamic,       bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::RangeBuffered,      bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::RangeFullyBuffered, bool, uint32_t, uint64_t >;
  template class CRSMatrix< crs_matrix::RangeCompressed,    bool, uint32_t, uint64_t >;

}  /* --- end of namespace diverg --- */
