# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25...3.30)

if (TARGET raptor::upgrade)
    return ()
endif ()

add_library ("raptor_upgrade" STATIC raptor_upgrade.cpp)
target_link_libraries ("raptor_upgrade" PUBLIC "raptor::interface")
add_library (raptor::upgrade ALIAS raptor_upgrade)
