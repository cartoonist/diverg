# Enabling `assert` by setting the build type to 'Debug' for test module
set(CMAKE_BUILD_TYPE "Debug")

# Prepare "Catch" library for other executables
if(NOT TARGET Catch2::Catch2)
  set(CATCH_SOURCE_DIR ext/catch2)
  execute_process(COMMAND git submodule update --init -- ${CATCH_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${CATCH_SOURCE_DIR})
endif()

# Defining the configure header file
configure_file(src/test_config.hpp.in ${PROJECT_BINARY_DIR}/test/src/test_config.hpp @ONLY)

# Default compile option
set(DIVERG_TESTS_DEFAULT_CXXOPS -g -Wall -Wpedantic -Wno-c++20-extensions)
# Default link option
set(DIVERG_TESTS_DEFAULT_LDOPS)

if(DIVERG_STRICT_ON_WARNS)
  list(APPEND DIVERG_TESTS_DEFAULT_CXXOPS -Werror)
endif(DIVERG_STRICT_ON_WARNS)

# Defining target diverg-tests
file(GLOB TEST_SOURCES "src/*.cpp")
add_executable(diverg-tests ${TEST_SOURCES})
target_compile_options(diverg-tests PRIVATE ${DIVERG_TESTS_DEFAULT_CXXOPS})
target_link_options(diverg-tests PRIVATE ${DIVERG_TESTS_DEFAULT_LDOPS})
target_include_directories(diverg-tests
  PRIVATE ${PROJECT_BINARY_DIR}/test/src
  # TODO: use diverg::libdiverg
  PRIVATE diverg::libdiverg_headers)
target_link_libraries(diverg-tests
  PRIVATE Catch2::Catch2
  # TODO: use diverg::libdiverg
  PRIVATE diverg::libdiverg_headers)
