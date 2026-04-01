# Recommanded CMake version 3.21
cmake_minimum_required(VERSION 3.16)

# Options
option(BUILD_TESTING "Build test programs" OFF)
option(BUILD_BENCHMARKING "Build benchmark program" OFF)
option(DIVERG_ETI "Build compiled library with ETI" OFF)
option(DIVERG_ENABLE_OPENMP "Enable OpenMP as host execution space" ON)
option(DIVERG_ENABLE_CUDA "Enable CUDA as device execution space" OFF)
option(DIVERG_STRICT_ON_WARNS "Pass '-Werror' flag to compiler" ON)
option(DIVERG_STATS "Enable statistics collection" OFF)
option(USE_BUNDLED_KOKKOS "Use bundled Kokkos" OFF)
option(USE_BUNDLED_KOKKOS_KERNELS "Use bundled KokkosKernels" OFF)
option(USE_BUNDLED_GUM "Use bundled gum library" OFF)
option(USE_BUNDLED_ALL "Use all bundled dependencies" OFF)

if(USE_BUNDLED_ALL)
  set(USE_BUNDLED_KOKKOS on)
  set(USE_BUNDLED_KOKKOS_KERNELS on)
  set(USE_BUNDLED_GUM on)
endif(USE_BUNDLED_ALL)

project(diverg VERSION 0.9.0 LANGUAGES CXX)

# Include external modules
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Finding dependencies.
if(DIVERG_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
endif(DIVERG_ENABLE_OPENMP)
find_package(Kokkos QUIET)
find_package(KokkosKernels QUIET)
find_package(gum 2.0.2 QUIET)
# Handle bundled dependencies
include(${PROJECT_SOURCE_DIR}/ext/external.cmake)

# Check if installed Kokkos has
if(DIVERG_ENABLE_OPENMP AND NOT Kokkos_ENABLE_OPENMP)
  message(FATAL "OpenMP is specified as target host execution space,"
    " but it is not enabled in Kokkos found by the build system")
endif()

if(DIVERG_ENABLE_CUDA AND NOT Kokkos_ENABLE_CUDA)
  message(FATAL "CUDA is specified as target device execution space,"
    " but it is not enabled in Kokkos found by the build system")
endif()

# NOTE: (temporarily) go easy on warnings when building with nvcc or for Release
if((CMAKE_BUILD_TYPE MATCHES Release) OR DIVERG_ENABLE_CUDA)
  set(DIVERG_STRICT_ON_WARNS off)
endif()

# Creating the INTERFACE library for header-only library
add_library(libdiverg_headers INTERFACE)
# Defining include directories
target_include_directories(libdiverg_headers
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>;$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE $<TARGET_PROPERTY:gum::gum,INTERFACE_INCLUDE_DIRECTORIES>
  INTERFACE $<TARGET_PROPERTY:Kokkos::kokkos,INTERFACE_INCLUDE_DIRECTORIES>
  INTERFACE $<TARGET_PROPERTY:Kokkos::kokkoskernels,INTERFACE_INCLUDE_DIRECTORIES>)
# Defining link libraries
target_link_libraries(libdiverg_headers
  INTERFACE $<BUILD_INTERFACE:gum::gum>;$<INSTALL_INTERFACE:gum::gum>
  INTERFACE $<BUILD_INTERFACE:Kokkos::kokkos>;$<INSTALL_INTERFACE:Kokkos::kokkos>
  INTERFACE $<BUILD_INTERFACE:Kokkos::kokkoskernels>;$<INSTALL_INTERFACE:Kokkos::kokkoskernels>)
# Use C++14
target_compile_features(libdiverg_headers INTERFACE cxx_std_14)

# Define `DIVERG_STATS`
if(DIVERG_STATS)
  target_compile_definitions(libdiverg_headers INTERFACE DIVERG_STATS)
endif(DIVERG_STATS)

# Generating the configure header file
configure_file(include/diverg/config.hpp.in ${PROJECT_SOURCE_DIR}/include/diverg/config.hpp @ONLY)

# Setting header files
file(GLOB HEADER_FILES RELATIVE "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/include/diverg/*.hpp")
# Keep ETI headers separate: they install into diverg/eti/
file(GLOB ETI_HEADER_FILES RELATIVE "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/include/diverg/eti/*.hpp")
list(TRANSFORM HEADER_FILES PREPEND "${PROJECT_SOURCE_DIR}/include/" OUTPUT_VARIABLE BUILD_HEADER_FILES)
list(TRANSFORM ETI_HEADER_FILES PREPEND "${PROJECT_SOURCE_DIR}/include/" OUTPUT_VARIABLE BUILD_ETI_HEADER_FILES)
list(TRANSFORM HEADER_FILES PREPEND "${CMAKE_INSTALL_FULL_INCLUDEDIR}/" OUTPUT_VARIABLE INSTALL_HEADER_FILES)
list(TRANSFORM ETI_HEADER_FILES PREPEND "${CMAKE_INSTALL_FULL_INCLUDEDIR}/" OUTPUT_VARIABLE INSTALL_ETI_HEADER_FILES)
# Defining target sources (ETI headers included for IDE visibility)
target_sources(libdiverg_headers
  INTERFACE "$<BUILD_INTERFACE:${BUILD_HEADER_FILES};${BUILD_ETI_HEADER_FILES}>;$<INSTALL_INTERFACE:${INSTALL_HEADER_FILES};${INSTALL_ETI_HEADER_FILES}>")

# Defining `diverg::libdiverg_headers` alias
add_library(diverg::libdiverg_headers ALIAS libdiverg_headers)

# Creating the compiled library for ETI mode
if(DIVERG_ETI)
  set(DIVERG_ETI_SOURCES
    src/eti/hbitvector.cpp
    src/eti/crs_matrix.cpp
    src/eti/range_sparse.cpp)
  if(DIVERG_ENABLE_CUDA)
    list(APPEND DIVERG_ETI_SOURCES
      src/eti/hbitvector_cuda.cpp
      src/eti/range_sparse_cuda.cpp)
    # ETI sources include Kokkos headers with KOKKOS_FUNCTION decorators
    # which expand to `__device__ __host__` which requires a CUDA compiler.
    # Enable the CUDA language in CMake if not already active (find_package(Kokkos)
    # may not do this automatically for downstream projects).
    enable_language(CUDA)
    set_source_files_properties(
      src/eti/hbitvector.cpp
      src/eti/hbitvector_cuda.cpp
      src/eti/range_sparse.cpp
      src/eti/range_sparse_cuda.cpp
      PROPERTIES LANGUAGE CUDA)
  endif()
  add_library(libdiverg ${DIVERG_ETI_SOURCES})
  # Set the output name for the library target
  set_target_properties(libdiverg PROPERTIES OUTPUT_NAME "diverg")
  if(DIVERG_ENABLE_CUDA)
    # Device symbols from explicit instantiations must be linkable across TUs.
    set_target_properties(libdiverg PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  endif()
  # Inherit all include paths and link deps from the INTERFACE header-only target
  target_link_libraries(libdiverg PUBLIC libdiverg_headers)
  # Advertise ETI mode to users so they see the extern template declarations
  target_compile_definitions(libdiverg PUBLIC DIVERG_ETI_ENABLED)
  add_library(diverg::libdiverg ALIAS libdiverg)
  add_library(diverg::diverg ALIAS libdiverg)
else()
  add_library(diverg::diverg ALIAS libdiverg_headers)
endif()

# Install headers
install(FILES ${BUILD_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/diverg)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/diverg/eti
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/diverg)
# Export targets
if(DIVERG_ETI)
  install(TARGETS libdiverg_headers libdiverg EXPORT diverg-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
else()
  install(TARGETS libdiverg_headers EXPORT diverg-targets)
endif()
install(EXPORT diverg-targets NAMESPACE diverg:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diverg)
# Creating the package configuration file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/diverg-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/diverg-config.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diverg)
# Generating the version file for the configuration file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/diverg-config-version.cmake"
  VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}"
  COMPATIBILITY AnyNewerVersion)
# Install generated configuration files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/diverg-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/diverg-config-version.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/diverg)

# Adding benchmark submodule
if(BUILD_BENCHMARKING)
  add_subdirectory(benchmark)
endif(BUILD_BENCHMARKING)

# Adding test submodule
if(BUILD_TESTING)
  add_subdirectory(test)
endif(BUILD_TESTING)

# Registering tests.
include(CTest)
add_test(NAME TestHBitVector COMMAND diverg-tests "[hbitvector]")
add_test(NAME TestRangeSparse COMMAND diverg-tests "[range_sparse]")
add_test(NAME TestCRSMatrix COMMAND diverg-tests "[crsmatrix]")
