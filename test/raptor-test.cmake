# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

# ----------------------------------------------------------------------------
# Short-circuit if tests are already configured
# ----------------------------------------------------------------------------

if (TARGET raptor::test)
    return ()
endif ()

message (STATUS "${ColourBold}Configuring tests${ColourReset}")

# ----------------------------------------------------------------------------
# Add Raptor
# ----------------------------------------------------------------------------

get_filename_component (Raptor_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
include ("${Raptor_SOURCE_DIR}/cmake/version.cmake")
include ("${Raptor_SOURCE_DIR}/cmake/configuration.cmake")
add_subdirectory ("${Raptor_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/raptor")
target_compile_options (raptor_interface INTERFACE "-pedantic" "-Wall" "-Wextra")

option (RAPTOR_WITH_WERROR "Report compiler warnings as errors." ON)

if (RAPTOR_WITH_WERROR)
    target_compile_options (raptor_interface INTERFACE "-Werror")
    message (STATUS "Building tests with -Werror.")
endif ()

# ----------------------------------------------------------------------------
# CPM
# ----------------------------------------------------------------------------

set (CPM_INDENT "CMake Package Manager CPM: ")
CPMUsePackageLock ("${CMAKE_CURRENT_LIST_DIR}/../cmake/package-lock.cmake")

# ----------------------------------------------------------------------------
# Paths to cmake modules
# ----------------------------------------------------------------------------

find_path (RAPTOR_TEST_CMAKE_MODULE_DIR
           NAMES declare_datasource.cmake
           HINTS "${CMAKE_CURRENT_LIST_DIR}/cmake/"
)
list (APPEND CMAKE_MODULE_PATH "${RAPTOR_TEST_CMAKE_MODULE_DIR}")

# ----------------------------------------------------------------------------
# Interface targets for the different test modules
# ----------------------------------------------------------------------------

enable_testing ()

# ----------------------------------------------------------------------------
# raptor::test
# ----------------------------------------------------------------------------

add_library (raptor_test INTERFACE)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # GCC12 and above: Disable warning about std::hardware_destructive_interference_size not being ABI-stable.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12)
        target_compile_options (raptor_test INTERFACE "-Wno-interference-size")
    endif ()

    # Warn about failed return value optimization.
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 14)
        target_compile_options (raptor_interface INTERFACE "-Wnrvo")
    endif ()
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # std::views::join is experimental in LLVM 17
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18)
        target_compile_definitions (raptor_test INTERFACE "_LIBCPP_ENABLE_EXPERIMENTAL")
    endif ()
endif ()

target_link_libraries (raptor_test INTERFACE "raptor_lib")
target_include_directories (raptor_test INTERFACE "${CMAKE_CURRENT_LIST_DIR}/include")
add_library (raptor::test ALIAS raptor_test)

# ----------------------------------------------------------------------------
# raptor::test::performance
# ----------------------------------------------------------------------------

add_library (raptor_test_performance INTERFACE)
target_link_libraries (raptor_test_performance INTERFACE "raptor::test" "benchmark::benchmark_main")
add_library (raptor::test::performance ALIAS raptor_test_performance)

# ----------------------------------------------------------------------------
# raptor::test::unit
# ----------------------------------------------------------------------------

add_library (raptor_test_unit INTERFACE)
# GCC12 has some bogus warnings. They will not be fixed in googletest.
# https://github.com/google/googletest/issues/4232
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 12 AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 13)
        target_compile_options (raptor_test INTERFACE "-Wno-restrict")
    endif ()
endif ()
target_link_libraries (raptor_test_unit INTERFACE "raptor::test" "GTest::gtest_main")
add_library (raptor::test::unit ALIAS raptor_test_unit)

# ----------------------------------------------------------------------------
# raptor::test::header
# ----------------------------------------------------------------------------

add_library (raptor_test_header INTERFACE)
target_link_libraries (raptor_test_header INTERFACE "raptor::test::unit")
target_link_libraries (raptor_test_header INTERFACE "raptor::test::performance")
target_compile_options (raptor_test_header INTERFACE "-Wno-unused-function" "-Wno-unused-const-variable")
add_library (raptor::test::header ALIAS raptor_test_header)

# ----------------------------------------------------------------------------
# Commonly used macros for the different test modules in seqan3.
# ----------------------------------------------------------------------------

include (raptor_add_benchmark)
include (raptor_add_unit_test)
include (${CMAKE_CURRENT_LIST_DIR}/data/datasources.cmake)

# ----------------------------------------------------------------------------
# Set directories for test output files, input data and binaries.
# ----------------------------------------------------------------------------

file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)
add_definitions (-DOUTPUTDIR=\"${CMAKE_CURRENT_BINARY_DIR}/output/\")
add_definitions (-DDATADIR=\"${CMAKE_CURRENT_BINARY_DIR}/data/\")
add_definitions (-DBINDIR=\"${CMAKE_BINARY_DIR}/bin/\")
