/**
 *    @file  eti/range_sparse.hpp
 *   @brief  Extern template declarations for range sparse algorithm ETI.
 *
 *  This header is included automatically at the bottom of range_sparse.hpp
 *  unless DIVERG_RANGE_SPARSE_ETI_INST__ is defined (which the ETI source
 *  files define to suppress the externs before providing definitions).
 *
 *  ETI covers the high-level CRSMatrix-based algorithm overloads wrapped in
 *  RangeSpAdd, RangeSpGEMM, and RangePower class templates. The default
 *  partition for each accumulator is instantiated for grid::Auto,
 *  grid::Suggested, and grid::RunTime. grid::Fixed is excluded since it is
 *  impractical for ETI (infinite template parameter combinations).
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

#ifndef DIVERG_ETI_RANGE_SPARSE_HPP__
#define DIVERG_ETI_RANGE_SPARSE_HPP__

#include <cstdint>
#include <diverg/eti_macros.hpp>

/* === Helper macros (they are #undef'd at the end) === */

/* --- BTree config for a given grid and execution space --- */
#define _DIVERG_ETI_BTREE(GRID, EXEC) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, BTreeAccumulator, ThreadRangePolicyPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, BTreeAccumulator, ThreadRangePolicyPartition, EXEC > > )

/* --- HBitVector configs for all host L1 sizes (1024-–32768) --- */
#define _DIVERG_ETI_HBV_HOST(GRID, EXEC) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<32768>, TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<32768>, TeamSequentialPartition, EXEC > > )

/* --- HBitVector configs for GPU L1 sizes (1024-–16384, no 32768) --- */
#define _DIVERG_ETI_HBV_GPU(GRID, EXEC) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangeSpGEMM< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<1024>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<2048>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<4096>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<8192>,  TeamSequentialPartition, EXEC > > ) \
  DIVERG_ETI_EXTERN_CLASS( RangePower< _eti_rcrs_t, \
      SparseConfig< GRID, HBitVectorAccumulator<16384>, TeamSequentialPartition, EXEC > > )

/* --- All host-space configs for a given grid spec --- */
#define _DIVERG_ETI_HOST(GRID, EXEC) \
  _DIVERG_ETI_BTREE(GRID, EXEC) \
  _DIVERG_ETI_HBV_HOST(GRID, EXEC)

namespace diverg {

  /* === Primary matrix type for ETI === */
  using _eti_rcrs_t = CRSMatrix< crs_matrix::RangeDynamic, bool, uint32_t, uint64_t >;

  /* === Kokkos::Serial === */
#if defined(KOKKOS_ENABLE_SERIAL)
  DIVERG_ETI_EXTERN_CLASS( RangeSpAdd< _eti_rcrs_t, Kokkos::Serial > )
  _DIVERG_ETI_HOST(grid::Auto,      Kokkos::Serial)
  _DIVERG_ETI_HOST(grid::Suggested,  Kokkos::Serial)
  _DIVERG_ETI_HOST(grid::RunTime,    Kokkos::Serial)
#endif

  /* === Kokkos::OpenMP === */
#if defined(KOKKOS_ENABLE_OPENMP)
  DIVERG_ETI_EXTERN_CLASS( RangeSpAdd< _eti_rcrs_t, Kokkos::OpenMP > )
  _DIVERG_ETI_HOST(grid::Auto,      Kokkos::OpenMP)
  _DIVERG_ETI_HOST(grid::Suggested,  Kokkos::OpenMP)
  _DIVERG_ETI_HOST(grid::RunTime,    Kokkos::OpenMP)
#endif

  /* === Kokkos::Cuda  (no BTree — host-only; max L1 = 16384 for GPU) === */
#if defined(KOKKOS_ENABLE_CUDA)
  DIVERG_ETI_EXTERN_CLASS( RangeSpAdd< _eti_rcrs_t, Kokkos::Cuda > )
  _DIVERG_ETI_HBV_GPU(grid::Auto,      Kokkos::Cuda)
  _DIVERG_ETI_HBV_GPU(grid::Suggested,  Kokkos::Cuda)
  _DIVERG_ETI_HBV_GPU(grid::RunTime,    Kokkos::Cuda)
#endif

}  /* --- end of namespace diverg --- */

#undef _DIVERG_ETI_BTREE
#undef _DIVERG_ETI_HBV_HOST
#undef _DIVERG_ETI_HBV_GPU
#undef _DIVERG_ETI_HOST

#endif  /* --- #ifndef DIVERG_ETI_RANGE_SPARSE_HPP__ --- */
