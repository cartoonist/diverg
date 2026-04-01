/**
 *    @file  src/eti/range_sparse.cpp
 *   @brief  Explicit template instantiations for range sparse algorithm
 *           wrapper classes (host execution spaces).
 *
 *  This file covers host version of ETI definitions (Kokkos::Serial and
 *  Kokkos::OpenMP).
 *  
 *  For device version (Kokkos::Cuda) see `range_sparse_cuda.cpp`.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  @internal
 *       Created:  Wed Apr 02, 2026  00:00
 *  Organization:  Universitat Bielefeld
 *     Copyright:  Copyright (c) 2026, Ali Ghaffaari
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

/* NOTE: Define this guard BEFORE including `range_sparse.hpp` so that it skips
    the extern declarations (this source file provides the definitions). */
#define DIVERG_RANGE_SPARSE_ETI_INST__

#include <diverg/range_sparse.hpp>

/* === Helper macros === */

#define _DIVERG_INST_BTREE(GRID, EXEC) \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, BTreeAccumulator, ThreadRangePolicyPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, BTreeAccumulator, ThreadRangePolicyPartition, EXEC > >;

#define _DIVERG_INST_HBV_HOST(GRID, EXEC) \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > >; \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > >; \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > >; \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > >; \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > >; \
  template class RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<32768>, TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<32768>, TeamSequentialPartition, EXEC > >;

#define _DIVERG_INST_HOST(GRID, EXEC) \
  _DIVERG_INST_BTREE(GRID, EXEC) \
  _DIVERG_INST_HBV_HOST(GRID, EXEC)

namespace diverg {

  using _eti_rcrs_t = CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t >;

#if defined(KOKKOS_ENABLE_SERIAL)
  template class RangeSpAdd< _eti_rcrs_t, Kokkos::Serial >;
  _DIVERG_INST_HOST(grid::Auto,      Kokkos::Serial)
  _DIVERG_INST_HOST(grid::Suggested,  Kokkos::Serial)
  _DIVERG_INST_HOST(grid::RunTime,    Kokkos::Serial)
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
  template class RangeSpAdd< _eti_rcrs_t, Kokkos::OpenMP >;
  _DIVERG_INST_HOST(grid::Auto,      Kokkos::OpenMP)
  _DIVERG_INST_HOST(grid::Suggested,  Kokkos::OpenMP)
  _DIVERG_INST_HOST(grid::RunTime,    Kokkos::OpenMP)
#endif

}  /* --- end of namespace diverg --- */

#undef _DIVERG_INST_BTREE
#undef _DIVERG_INST_HBV_HOST
#undef _DIVERG_INST_HOST
