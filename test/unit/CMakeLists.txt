# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

project (raptor_test_unit LANGUAGES CXX C)

include (../raptor-test.cmake)
include (GoogleTest OPTIONAL)

CPMGetPackage (GTest)

add_subdirectory (api)
add_subdirectory (cli)
