# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.25)

if (TARGET raptor::update)
    return ()
endif ()

add_library ("raptor_update" STATIC raptor_update.cpp delete_user_bins.cpp dump_index.cpp insert_user_bin.cpp
insert/get_location.cpp
insert/compute_kmers.cpp
insert/required_technical_bins.cpp
insert/max_ibf_sizes.cpp
insert/insert_into_ibf.cpp
insert/insert_tb_and_parents.cpp
insert/update_bookkeeping.cpp
insert/max_elements.cpp
insert/find_ibf_size_splitting.cpp
insert/find_empty_bin_idx.cpp)
target_link_libraries ("raptor_update" PUBLIC "raptor::interface")
add_library (raptor::update ALIAS raptor_update)
