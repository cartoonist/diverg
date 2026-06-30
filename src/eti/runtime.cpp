/**
 *    @file  src/eti/runtime.cpp
 *   @brief  Implementation of the Kokkos-free runtime wrapper.
 *
 *  This translation unit wraps the Kokkos runtime calls behind the Kokkos-free
 *  interface declared in `runtime.hpp`. It includes Kokkos, so it is compiled
 *  by the device compiler when CUDA is enabled and lives in the compiled
 *  `libdiverg` library.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#include <diverg/runtime.hpp>

#include <Kokkos_Core.hpp>


namespace diverg {
  bool
  is_runtime_initialized( ) noexcept
  {
    return Kokkos::is_initialized();
  }

  void
  initialize_runtime( )
  {
    if ( !Kokkos::is_initialized() ) Kokkos::initialize();
  }

  void
  finalize_runtime( )
  {
    if ( Kokkos::is_initialized() ) Kokkos::finalize();
  }

  ScopedRuntime::ScopedRuntime( bool finalize_on_destroy )
    : m_finalize( finalize_on_destroy )
  {
    if ( !Kokkos::is_initialized() ) Kokkos::initialize();
  }

  ScopedRuntime::~ScopedRuntime( )
  {
    if ( this->m_finalize && Kokkos::is_initialized() ) Kokkos::finalize();
  }
}  /* --- end of namespace diverg --- */
