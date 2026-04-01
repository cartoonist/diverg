/**
 *    @file  src/eti/range_sparse_cuda.cpp
 *   @brief  Explicit template instantiations for range sparse algorithm
 *           wrapper classes (Kokkos::Cuda execution space).
 *
 *  This file must be compiled by the CUDA compiler (nvcc).
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

#if defined(KOKKOS_ENABLE_CUDA)

/* === Helper macros === */

#define _DIVERG_INST_HBV_GPU(GRID, EXEC) \
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
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > >; \
  template class RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > >;

namespace diverg {

  using _eti_rcrs_t = CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t >;

  template class RangeSpAdd< _eti_rcrs_t, Kokkos::Cuda >;
  _DIVERG_INST_HBV_GPU(grid::Auto,      Kokkos::Cuda)
  _DIVERG_INST_HBV_GPU(grid::Suggested,  Kokkos::Cuda)
  _DIVERG_INST_HBV_GPU(grid::RunTime,    Kokkos::Cuda)

}  /* --- end of namespace diverg --- */

#undef _DIVERG_INST_HBV_GPU

#endif  /* KOKKOS_ENABLE_CUDA */
