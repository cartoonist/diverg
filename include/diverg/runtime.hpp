/**
 *    @file  runtime.hpp
 *   @brief  Kokkos-free runtime (execution space) management wrapper.
 *
 *  This header lets downstream code initialise and finalise the Kokkos runtime
 *  without including any Kokkos header. In an ETI-build, the implementation
 *  lives in the compiled `libdiverg` library (built by the device compiler when
 *  CUDA is enabled). This makes it possible for a host translation unit to be
 *  compiled without any dependency on CUDA toolchain.
 *
 *  @author  Ali Ghaffaari (\@cartoonist), <ali.ghaffaari@uni-bielefeld.de>
 *
 *  This source code is released under the terms of the MIT License.
 *  See LICENSE file for more information.
 */

#ifndef DIVERG_RUNTIME_HPP__
#define DIVERG_RUNTIME_HPP__


namespace diverg {
  /**
   *  @brief  Whether the runtime has been initialised.
   */
  bool is_runtime_initialized( ) noexcept;

  /**
   *  @brief  Initialise the runtime with default settings (if not initialised).
   */
  void initialize_runtime( );

  /**
   *  @brief  Finalise the runtime (if it is initialised).
   */
  void finalize_runtime( );

  /**
   *  @brief  RAII guard initialising and optionally finalises the runtime.
   *
   *  The main purpose of this class is to manage the Kokkos runtime.
   */
  class ScopedRuntime {
  public:
    explicit ScopedRuntime( bool finalize_on_destroy = true );
    ~ScopedRuntime();

    ScopedRuntime( ScopedRuntime const& ) = delete;
    ScopedRuntime& operator=( ScopedRuntime const& ) = delete;
    ScopedRuntime( ScopedRuntime&& ) = delete;
    ScopedRuntime& operator=( ScopedRuntime&& ) = delete;

  private:
    bool m_finalize;
  };
}  // namespace diverg

#endif  /* --- #ifndef DIVERG_RUNTIME_HPP__ --- */
