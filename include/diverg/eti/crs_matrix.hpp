/**
 *    @file  eti/crs_matrix.hpp
 *   @brief  Extern template declarations for CRSMatrix ETI.
 *
 *  This header is included automatically at the bottom of crs_matrix.hpp
 *  unless DIVERG_CRSMATRIX_ETI_INST__ is defined (which the ETI source
 *  file defines to suppress the externs before providing definitions).
 *
 *  CRSMatrix does not depend on execution space; ordinal=uint32_t,
 *  size=uint64_t covers the primary use case. Extend as needed.
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

#ifndef DIVERG_ETI_CRSMATRIX_HPP__
#define DIVERG_ETI_CRSMATRIX_HPP__

#include <cstdint>
#include <diverg/eti_macros.hpp>

namespace diverg {

  /* === Basic group === */
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::Dynamic,       bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::Buffered,      bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::FullyBuffered, bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::Compressed,    bool, uint32_t, uint64_t> )

  /* === Range group === */
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::RangeDynamic,       bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::RangeBuffered,      bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::RangeFullyBuffered, bool, uint32_t, uint64_t> )
  DIVERG_ETI_EXTERN_CLASS( CRSMatrix<crs_matrix::RangeCompressed,    bool, uint32_t, uint64_t> )

}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_ETI_CRSMATRIX_HPP__ --- */
