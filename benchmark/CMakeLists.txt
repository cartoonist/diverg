# Include external modules
include(GNUInstallDirs)

# Setting the build type if it is not already set
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Default compile option
set(DIVERG_BENCHMARK_DEFAULT_CXXOPS -g -Wall -Wpedantic)
# Default link option
set(DIVERG_BENCHMARK_DEFAULT_LDOPS)

if(DIVERG_STRICT_ON_WARNS)
  list(APPEND DIVERG_BENCHMARK_DEFAULT_CXXOPS -Werror)
endif(DIVERG_STRICT_ON_WARNS)

if(DIVERG_ENABLE_CUDA)
  list(APPEND DIVERG_BENCHMARK_DEFAULT_CXXOPS -Wno-deprecated-gpu-targets)
endif(DIVERG_ENABLE_CUDA)

# Defining target 'rcrs_benchmark': Range CRS matrix operations benchmark
set(RCRS_BENCH_SOURCES "src/rcrs_benchmark.cpp")
add_executable(rcrs_benchmark ${RCRS_BENCH_SOURCES})
if(DIVERG_ENABLE_CUDA AND DIVERG_ETI)
  # ETI builds do not use nvcc_wrapper compile launcher to build the targets;
  # i.e. `Kokkos_ENABLE_COMPILE_AS_CMAKE_LANGUAGE` is on. So, every TU including
  # Kokkos must be compiled in CUDA language (nvcc) and built as relocatable
  # device code to match Kokkos/desul.
  #
  # Without ETI, Kokkos uses the nvcc_wrapper launcher and host TUs are already
  # handled through nvcc, so they should not be re-tagged.
  set_source_files_properties(${RCRS_BENCH_SOURCES} PROPERTIES LANGUAGE CUDA)
  set_target_properties(rcrs_benchmark PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
target_compile_options(rcrs_benchmark PRIVATE ${DIVERG_BENCHMARK_DEFAULT_CXXOPS})
target_link_options(rcrs_benchmark PRIVATE ${DIVERG_BENCHMARK_DEFAULT_LDOPS})
target_include_directories(rcrs_benchmark
  PRIVATE diverg::diverg)
target_link_libraries(rcrs_benchmark
  PRIVATE diverg::diverg)
# Install targets
install(TARGETS rcrs_benchmark DESTINATION ${CMAKE_INSTALL_BINDIR})
