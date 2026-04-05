DiVerG
======
Scalable Distance Index for Validation of Paired-End Alignments in Sequence Graphs
----------------------------------------------------------------------------------

This is an implementation of DiVerG method introduced in:

> Ali Ghaffaari, Alexander Schönhuth, and Tobias Marschall.<br>
> *DiVerG: Scalable Distance Index for Validation of Paired-End Alignments in Sequence Graphs.*
> In 25th International Conference on Algorithms for Bioinformatics (WABI 2025).
> Leibniz International Proceedings in Informatics (LIPIcs), Volume 344, pp. 10:1-10:24,
>
> [DOI: 10.4230/LIPIcs.WABI.2025.10](https://doi.org/10.4230/LIPIcs.WABI.2025.10)

## Table of Contents

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Installation and Build](#installation-and-build)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Testing](#testing)
- [Benchmarking](#benchmarking)
- [License](#license)

## Introduction

DiVerG is an indexing scheme based on [PairG](https://github.com/ParBLiSS/PairG)
that efficiently determines whether there is a path between any two loci in a
sequence graph that falls within a particular, statistically well-motivated
range $[d_1, d_2]$. This problem is motivated by the need in paired-end read
mapping to sequence graph to verify distance constraints between candidate
alignments.

DiVerG introduces a compact data structure for representing Boolean sparse matrices,
as well as a fast and scalable algorithm for computing matrix-matrix
multiplication and addition using the compressed representation on CUDA and
OpenMP backends.

It overcomes the limitations of PairG by exploiting the
extensive potential for improvements in terms of scalability and space
efficiency. As a consequence, it can process substantially larger datasets, such
as whole human genomes, which are unmanageable by PairG.
DiVerG offers faster index construction time and consistently faster query time
with gains proportional to the size of the underlying compact data structure.

The constructed index can be used for efficient distance querying between
two loci for constraint validation.

## Dependencies

DiVerG is a C++ header-only library which can also be compiled as a static
library using [Explicit Template Instantiation (ETI)](#explicit-template-instantiation-eti)
to reduce downstream compile times and/or decouple downstream translation units
from CUDA compiler requirements. It requires:

- **C++17** compatible compiler
- **CMake** 3.21 or later (3.16 minimum, 3.21 recommended)
- **[gum](https://github.com/cartoonist/gum)** (>= 2.0.2) -- sequence graph library
- **[Kokkos](https://github.com/kokkos/kokkos)** -- performance portability framework

Optional dependencies:

- **[KokkosKernels](https://github.com/kokkos/kokkos-kernels)** -- Required by
  utility functions in `range_sparse_utils.hpp` header.

> [!NOTE]
> The test and benchmark modules require the utilities in
> `range_sparse_utils.hpp` and enable CMake option `DIVERG_ENABLE_UTILS`
> automatically which is by default OFF and the header file is excluded in
> builds without these modules.

- **CUDA** - For GPU acceleration (requires CUDA-enabled Kokkos)
- **OpenMP** - For CPU parallelization (enabled by default)

The build system can fetch and build dependencies automatically by passing
`-DUSE_BUNDLED_<dependency>=ON` or `-DUSE_BUNDLED_ALL=ON` during CMake
configuration. Individual bundled dependency options:

- `USE_BUNDLED_GUM`
- `USE_BUNDLED_KOKKOS`
- `USE_BUNDLED_KOKKOS_KERNELS`

By turning these options OFF, the build system looks for installed versions of
dependencies and does not attempt to install third-party libraries in install
prefix (`CMAKE_INSTALL_PREFIX`).

## Installation and Build

DiVerG uses CMake with preset configurations for different build types.

### Quick Start

```bash
# Clone the repository
git clone https://github.com/cartoonist/diverg.git
cd diverg

# Release build (header-only) with OpenMP (fetches all dependencies)
cmake --preset default-release-omp
cmake --build build/Release
# OR
# Release build (compiled library) with CUDA + OpenMP
cmake --preset default-release-eti-cuda
cmake --build build/ReleaseCuda
```

### Build Presets

Presets are named `[default|ninja]-<TYPE>[-eti][-all]-<BACKEND>`. Each preset is
available with both GNU Make (prefix `default-`) and Ninja (prefix `ninja-`) as
generators.

- **TYPE**: `debug`, `release`, or `benchmark` (release with `DIVERG_STATS`)
- **eti**: Presets containing `-eti-` build the compiled ETI library
- **all**: Presets containing `-all-` include tests and benchmarks
- **BACKEND**: `omp` (OpenMP) or `cuda` (CUDA + OpenMP)

For example, `ninja-release-eti-cuda` configures a release build with the
compiled library using Ninja build system.

### Build Options

| Option                   | Default | Description                                       |
|--------------------------|---------|---------------------------------------------------|
| `BUILD_TESTING`          | OFF     | Build test programs (requires KokkosKernels)      |
| `BUILD_BENCHMARKING`     | OFF     | Build benchmark program (requires KokkosKernels)  |
| `DIVERG_ETI`             | OFF     | Build compiled library with ETI                   |
| `DIVERG_ENABLE_OPENMP`   | ON      | Enable OpenMP execution space                     |
| `DIVERG_ENABLE_CUDA`     | OFF     | Enable CUDA execution space                       |
| `DIVERG_ENABLE_UTILS`    | OFF     | Enable utility functions (requires KokkosKernels) |
| `DIVERG_STATS`           | OFF     | Enable statistics collection                      |
| `DIVERG_STRICT_ON_WARNS` | OFF     | Treat warnings as errors (`-Werror`)              |
| `USE_BUNDLED_ALL`        | OFF     | Fetch and build all dependencies                  |

## Usage

### As a Header-Only Library

Include the main header in your C++ project:

```cpp
#include <diverg/dindex.hpp>
```

### CMake Integration

After installing DiVerG, use `find_package` in your `CMakeLists.txt`:

```cmake
find_package(diverg REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE diverg::diverg)
```

The `diverg::diverg` target automatically resolves to the compiled ETI library
if it was installed, otherwise it falls back to the header-only target
(`diverg::libdiverg_headers`).

### Explicit Template Instantiation (ETI)

By default, DiVerG is header-only: all template definitions are included and
instantiated in every translation unit that uses them. This is simple but can
increase compile times for downstream projects and most importantly, when using
CUDA, introduces device functions to your source code which may create variety
of (unexpected) issues when compiling your code.

Enabling `DIVERG_ETI=ON` builds a compiled library (`libdiverg`) that contains
pre-instantiated templates for the common type combinations.

ETI presets are provided for both OpenMP and CUDA backends:

```bash
# Install with ETI
cmake --preset ninja-release-eti-omp
cmake --build build/ReleaseEti
cmake --install build/ReleaseEti --prefix /path/to/install
```

## Project Structure

```
include/diverg/       Header-only library implementation
  ├── dindex.hpp          Distance index methods (main entry point)
  ├── range_sparse.hpp    Range sparse matrix algorithms and data structures
  ├── hbitvector.hpp      Hierarchical bit vector
  ├── crs_matrix.hpp      CRS and rCRS matrix definitions and specialisations
  ├── eti_macros.hpp      ETI macro definitions
  └── eti/                Explicit template instantiation headers
src/eti/              ETI source files (compiled when DIVERG_ETI=ON)
test/                 Test suite (Catch2) and test data
benchmark/            Performance benchmarks
cmake/                CMake package configuration
ext/                  External dependency build scripts
CMakePresets.json     CMake preset configurations
```

## Testing

Tests use [Catch2](https://github.com/catchorg/Catch2) and are built
automatically with the debug presets (`BUILD_TESTING=ON`).

```bash
# Build and run all tests
cmake --preset default-debug-all-omp
cmake --build build/Debug  # with -j if you are using GNU Make
cmake --build build/Debug --target test

# Run specific test categories
./build/Debug/test/diverg-tests "[hbitvector]"
./build/Debug/test/diverg-tests "[range_sparse]"
./build/Debug/test/diverg-tests "[crsmatrix]"

# Using CTest
cd build/Debug && ctest
```

## Benchmarking

```bash
# Configure and build
cmake --preset ninja-benchmark-all-omp
cmake --build build/Benchmark

# Run range sparse matrix benchmark
./build/Benchmark/benchmark/rcrs_benchmark [ARGS]
```

The benchmark preset enables `DIVERG_STATS` for collecting performance
statistics during execution.

## License

This project is licensed under the MIT License -- see the [LICENSE](LICENSE) file for details.
