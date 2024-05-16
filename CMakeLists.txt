# Recommanded CMake version 3.21
cmake_minimum_required(VERSION 3.16)

# Options
option(BUILD_TESTING "Build test programs" OFF)
option(BUILD_BENCHMARKING "Build benchmark program" OFF)
# TODO: handle BUILD_HEADERS_ONLY option
option(BUILD_HEADERS_ONLY "Configure as header-only library" OFF)
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
list(TRANSFORM HEADER_FILES PREPEND "${PROJECT_SOURCE_DIR}/include/" OUTPUT_VARIABLE BUILD_HEADER_FILES)
list(TRANSFORM HEADER_FILES PREPEND "${CMAKE_INSTALL_FULL_INCLUDEDIR}/" OUTPUT_VARIABLE INSTALL_HEADER_FILES)
# Defining target sources
target_sources(libdiverg_headers
  INTERFACE "$<BUILD_INTERFACE:${BUILD_HEADER_FILES}>;$<INSTALL_INTERFACE:${INSTALL_HEADER_FILES}>")

# Defining `diverg::libdiverg_headers` alias
add_library(diverg::libdiverg_headers ALIAS libdiverg_headers)

# Install targets for libdiverg_headers
install(FILES ${BUILD_HEADER_FILES} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/diverg)
# Exporting targets
install(TARGETS libdiverg_headers EXPORT diverg-targets)
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
