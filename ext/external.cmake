# When `kokkos` is not found
if(NOT TARGET Kokkos::kokkos)
  if(NOT USE_BUNDLED_KOKKOS)
    message(FATAL_ERROR "Kokkos library not found. "
      "Pass in `-DUSE_BUNDLED_KOKKOS=on` when running cmake to use the bundled version. "
      "It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled Kokkos library")
  set(Kokkos_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kokkos)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${Kokkos_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  if(DIVERG_ENABLE_OPENMP)
    set(Kokkos_ENABLE_OPENMP on CACHE BOOL "Enable Kokkos OpenMP space")
  endif()
  if(DIVERG_ENABLE_CUDA)
    set(Kokkos_ENABLE_CUDA on CACHE BOOL "Enable Kokkos CUDA space")
    if(DIVERG_ETI)
      # The ETI library is built with relocatable device code (RDC)
      # (CUDA_SEPARABLE_COMPILATION) so the device symbols from ETIs can be
      # linked across translation units. Kokkos (and its dependencies; e.g.
      # desul atomics) must be configured for RDC to match.
      #
      # These following settings are required for a correct DiVerG CUDA build,
      # so they are FORCEd to override Kokkos defaults.
      set(Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE on CACHE BOOL
          "Enable Kokkos CUDA relocatable device code" FORCE)
      # Build Kokkos against CMake's native CUDA language instead of installing
      # the global nvcc_wrapper launcher. With nvcc_wrapper, EVERY target
      # including Kokkos is routed through nvcc which may breaks downstream host
      # code that cannot be compiled by nvcc. In native-language mode only the
      # CUDA-language ETI sources are compiled with nvcc, host code keeps using
      # the host compiler, and Kokkos propagates its required nvcc options to
      # CUDA-language targets.
      set(Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE on CACHE BOOL
          "Use the native CMake CUDA language instead of the nvcc_wrapper launcher" FORCE)
    endif()
  endif()
  add_subdirectory(${Kokkos_SOURCE_DIR})
endif()

# When `kokkos-kernels` is not found (optional: only when DIVERG_ENABLE_UTILS)
if(DIVERG_ENABLE_UTILS AND NOT TARGET Kokkos::kokkoskernels)
  if(NOT USE_BUNDLED_KOKKOS_KERNELS)
    message(FATAL_ERROR "KokkosKernels library not found. "
      "Pass in `-DUSE_BUNDLED_KOKKOS_KERNELS=on` when running cmake to use the bundled "
      "version. It will be installed alongside the library. "
      "Alternatively, pass `-DDIVERG_ENABLE_UTILS=off` to build without it.")
  endif()
  message(STATUS "Using bundled kokkos-kernels library")
  set(KokkosKernels_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kokkos-kernels)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KokkosKernels_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  if(Kokkos_ENABLE_CUDA AND Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
    # When embedded (KOKKOSKERNELS_HAS_PARENT), KokkosKernels skips its
    # standalone `find_package(Kokkos)` branch and therefore neither receives
    # the `Kokkos_COMPILE_LANGUAGE` variable nor enables the compile language
    # itself when Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE=TRUE.
    # Without these, KokkosKernels' device sources are considered as host C++
    # (which fail to be compiled). Provide both here (mirroring what
    # KokkosKernels' standalone branch does) before adding it.
    set(Kokkos_COMPILE_LANGUAGE CUDA)
    enable_language(CUDA)
    # In ETI builds, DiVerG forces relocatable device code (RDC) on, so desul
    # (bundled in Kokkos) is configured for separable compilation. A non-RDC
    # CUDA TU then trips its consistency `#error`. The embedded KokkosKernels
    # targets do not inherit RDC, so compiling the CUDA targets in its subtree
    # with `-rdc=true` is necessary.
    set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
  endif()
  add_subdirectory(${KokkosKernels_SOURCE_DIR})
endif()

# When `gum` is not found
if(NOT TARGET gum::gum)
  if(NOT USE_BUNDLED_GUM)
    message(FATAL_ERROR "gum library not found. "
      "Pass in `-DUSE_BUNDLED_GUM=on` when running cmake to use the bundled version. "
      "It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled gum library")
  set(GUM_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/gum)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${GUM_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set(BUILD_TESTING_SAVED "${BUILD_TESTING}")
  set(BUILD_TESTING OFF)
  # DiVerG only supports GFA.
  # NOTE: Downstream projects (e.g. PSI) accessing gum library can request them
  # by setting these before adding DiVerG as a subdirectory.
  if(NOT DEFINED GUM_WITH_VG)
    set(GUM_WITH_VG OFF)
  endif()
  if(NOT DEFINED GUM_WITH_VGIO)
    set(GUM_WITH_VGIO OFF)
  endif()
  if(NOT DEFINED GUM_WITH_HG)
    set(GUM_WITH_HG OFF)
  endif()
  if(NOT DEFINED GUM_WITH_BDSG)
    set(GUM_WITH_BDSG OFF)
  endif()
  set(GUM_USE_VCPKG OFF)
  add_subdirectory(${GUM_SOURCE_DIR})
  set(BUILD_TESTING "${BUILD_TESTING_SAVED}")
endif()
