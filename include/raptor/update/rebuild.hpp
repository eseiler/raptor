// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/argument_parsing/build_arguments.hpp>
#include <raptor/build/hibf/build_data.hpp>
#include <raptor/index.hpp>
#include <chopper/configuration.hpp>
#include <robin_hood.h>


namespace raptor
{
chopper::configuration layout_config(raptor_index<index_structure::hibf> & index,
                                      update_arguments const & arguments,
                                      std::string file_indicator="");

void call_layout(std::vector<std::tuple<size_t, std::string>> kmer_count_filenames, chopper::configuration & arguments);

std::vector<std::tuple<size_t, std::string>> get_kmer_counts(raptor_index<index_structure::hibf> const & index,
                                                             std::set<std::string> const & filenames);

void write_kmer_counts(std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
                       std::filesystem::path count_filename);

std::vector<uint64_t> split_mb(std::tuple<size_t,size_t> index_tuple,
           raptor_index<index_structure::hibf> & index,
           update_arguments const & update_arguments,
           int number_of_splits = 2);

std::vector<std::vector<std::tuple<size_t, std::string>>> split_filenames(
        std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
        int number_of_splits = 2);


void partial_rebuild(std::tuple<size_t,size_t> index_tuple,
                     raptor_index<index_structure::hibf> & index,
                     update_arguments const & arguments,
                     int number_of_splits = 2);


void full_rebuild(raptor_index<index_structure::hibf> & index,
                  update_arguments const & update_arguments);

void write_filenames(std::string bin_path, std::set<std::string> user_bin_filenames);

build_arguments build_config(raptor_index<index_structure::hibf> & index,
                             chopper::configuration layout_arguments);

template <seqan3::data_layout data_layout_mode>
robin_hood::unordered_flat_set<size_t> call_build(build_arguments & arguments,
                raptor_index<hierarchical_interleaved_bloom_filter<data_layout_mode>> & index,
                bool is_root = false);

template <typename T> void remove_indices(std::unordered_set<size_t> indices_to_remove, std::vector<T> & vector);

void remove_ibfs(raptor_index<index_structure::hibf> & index,
                 size_t ibf_idx);

void attach_subindex(raptor_index<index_structure::hibf> & index,
                   raptor_index<index_structure::hibf> & subindex,
                   std::tuple<size_t, size_t> index_tuple);

bool check_tmax_rebuild(raptor_index<index_structure::hibf> & index, size_t ibf_idx,
                        update_arguments const & update_arguments);

} // namespace raptor
