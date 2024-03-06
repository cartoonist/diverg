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
  endif()
  add_subdirectory(${Kokkos_SOURCE_DIR})
endif()

# When `kokkos-kernels` is not found
if(NOT TARGET Kokkos::kokkoskernels)
  if(NOT USE_BUNDLED_KOKKOS_KERNELS)
    message(FATAL_ERROR "KokkosKernels library not found. "
      "Pass in `-DUSE_BUNDLED_KOKKOS_KERNELS=on` when running cmake to use the bundled "
      "version. It will be installed alongside the library.")
  endif()
  message(STATUS "Using bundled kokkos-kernels library")
  set(KokkosKernels_SOURCE_DIR ${PROJECT_SOURCE_DIR}/ext/kokkos-kernels)
  execute_process(
    COMMAND git submodule update --init --recursive -- ${KokkosKernels_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
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
  set(GUM_WITH_VG OFF)
  set(GUM_WITH_VGIO OFF)
  set(GUM_WITH_HG OFF)
  set(GUM_WITH_BDSG OFF)
  set(GUM_USE_VCPKG OFF)
  add_subdirectory(${GUM_SOURCE_DIR})
  set(BUILD_TESTING "${BUILD_TESTING_SAVED}")
endif()
