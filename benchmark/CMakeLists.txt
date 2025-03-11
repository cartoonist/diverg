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
  list(APPEND DIVERG_TESTS_DEFAULT_CXXOPS -Wno-deprecated-gpu-targets)
endif(DIVERG_ENABLE_CUDA)

# Defining target 'rcrs_benchmark': Range CRS matrix operations benchmark
set(RCRS_BENCH_SOURCES "src/rcrs_benchmark.cpp")
add_executable(rcrs_benchmark ${RCRS_BENCH_SOURCES})
target_compile_options(rcrs_benchmark PRIVATE ${DIVERG_BENCHMARK_DEFAULT_CXXOPS})
target_link_options(rcrs_benchmark PRIVATE ${DIVERG_BENCHMARK_DEFAULT_LDOPS})
target_include_directories(rcrs_benchmark
  # TODO: Use `diverg::libdiverg`
  #PRIVATE diverg::libdiverg
  PRIVATE diverg::libdiverg_headers)
target_link_libraries(rcrs_benchmark
  # TODO: Use `diverg::libdiverg`
  #PRIVATE diverg::libdiverg
  PRIVATE diverg::libdiverg_headers)
# Install targets
install(TARGETS rcrs_benchmark DESTINATION ${CMAKE_INSTALL_BINDIR})
