/**
 *    @file  eti/hbitvector.hpp
 *   @brief  Extern template declarations for HBitVector ETI.
 *
 *  This header is included automatically at the bottom of hbitvector.hpp
 *  unless DIVERG_HBITVECTOR_ETI_INST__ is defined (which the ETI source
 *  file defines to suppress the externs before providing definitions).
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

#ifndef DIVERG_ETI_HBITVECTOR_HPP__
#define DIVERG_ETI_HBITVECTOR_HPP__

#include <diverg/eti_macros.hpp>

namespace diverg {

  /* === Kokkos::Serial === */
#if defined(KOKKOS_ENABLE_SERIAL)
  DIVERG_ETI_EXTERN_CLASS( HBitVector<1024,  Kokkos::Serial> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<2048,  Kokkos::Serial> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<4096,  Kokkos::Serial> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<8192,  Kokkos::Serial> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<16384, Kokkos::Serial> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<32768, Kokkos::Serial> )
#endif

  /* === Kokkos::OpenMP === */
#if defined(KOKKOS_ENABLE_OPENMP)
  DIVERG_ETI_EXTERN_CLASS( HBitVector<1024,  Kokkos::OpenMP> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<2048,  Kokkos::OpenMP> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<4096,  Kokkos::OpenMP> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<8192,  Kokkos::OpenMP> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<16384, Kokkos::OpenMP> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<32768, Kokkos::OpenMP> )
#endif

  /* === Kokkos::Cuda === */
  /* (GPU shared memory: max practical L1 = 16384 bits = 2 KB) */
#if defined(KOKKOS_ENABLE_CUDA)
  DIVERG_ETI_EXTERN_CLASS( HBitVector<1024,  Kokkos::Cuda> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<2048,  Kokkos::Cuda> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<4096,  Kokkos::Cuda> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<8192,  Kokkos::Cuda> )
  DIVERG_ETI_EXTERN_CLASS( HBitVector<16384, Kokkos::Cuda> )
#endif

}  /* --- end of namespace diverg --- */

#endif  /* --- #ifndef DIVERG_ETI_HBITVECTOR_HPP__ --- */
