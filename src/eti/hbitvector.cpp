/**
 *    @file  src/eti/hbitvector.cpp
 *   @brief  Explicit template instantiations for HBitVector (host execution spaces).
 *
 *  This file covers host version of ETI definitions (Kokkos::Serial and
 *  Kokkos::OpenMP).
 *  
 *  For device version (Kokkos::Cuda) see `hbitvector_cuda.cpp`.
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

/* NOTE: Define this guard BEFORE including `hbitvector.hpp` so that it skips
    the extern declarations (this source file provides the definitions). */
#define DIVERG_HBITVECTOR_ETI_INST__

#include <diverg/hbitvector.hpp>

namespace diverg {

#if defined(KOKKOS_ENABLE_SERIAL)
  template class HBitVector< 1024,  Kokkos::Serial >;
  template class HBitVector< 2048,  Kokkos::Serial >;
  template class HBitVector< 4096,  Kokkos::Serial >;
  template class HBitVector< 8192,  Kokkos::Serial >;
  template class HBitVector< 16384, Kokkos::Serial >;
  template class HBitVector< 32768, Kokkos::Serial >;
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  template class HBitVector< 1024,  Kokkos::OpenMP >;
  template class HBitVector< 2048,  Kokkos::OpenMP >;
  template class HBitVector< 4096,  Kokkos::OpenMP >;
  template class HBitVector< 8192,  Kokkos::OpenMP >;
  template class HBitVector< 16384, Kokkos::OpenMP >;
  template class HBitVector< 32768, Kokkos::OpenMP >;
#endif

}  /* --- end of namespace diverg --- */
