// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <hibf/build/insert_into_ibf.hpp>
#include <hibf/contrib/robin_hood.hpp>
#include <hibf/misc/divide_and_ceil.hpp>

#include <raptor/file_reader.hpp>
#include <raptor/update/insert_user_bin.hpp>

namespace raptor
{

struct bookkeeping_arguments
{
    size_t ibf_idx;
    size_t old_number_of_bins;
    size_t number_of_new_bins;
};

struct insert_location
{
    size_t ibf_idx;
    size_t bin_idx;
    size_t number_of_bins;
};

struct rebuild_location
{
    size_t ibf_idx;
    size_t bin_idx;
};

struct max_elements_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t bin_size{};
};

struct ibf_max
{
    size_t max_elements;
    size_t ibf_idx;

    constexpr auto operator<=>(ibf_max const & other) const = default;
};

struct required_technical_bins_parameters
{
    size_t bin_size{};
    size_t elements{};
    double fpr{};
    size_t hash_count{};
    size_t max_elements{};
};

struct ibf_location
{
    size_t ibf_idx;
    size_t max_elements;
};

robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::filesystem::path const & ub_file,
                                                       raptor_index<index_structure::hibf> const & index);

size_t
find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t const ibf_idx, size_t const number_of_bins);

ibf_location find_ibf_size_splitting(std::vector<ibf_max> const & max_ibf_sizes,
                                     size_t const kmer_count,
                                     raptor_index<index_structure::hibf> & index);

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index);

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index);

void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers,
                     insert_location const & insert_location,
                     raptor_index<index_structure::hibf> & index,
                     rebuild_location & rebuild_index_tuple);

rebuild_location insert_tb_and_parents(robin_hood::unordered_flat_set<size_t> const & kmers,
                                       insert_location insert_location,
                                       raptor_index<index_structure::hibf> & index);

size_t max_elements(max_elements_parameters const & params);

std::vector<ibf_max> max_ibf_sizes(raptor_index<index_structure::hibf> const & index);

size_t required_technical_bins(required_technical_bins_parameters const & params);

void update_bookkeeping(bookkeeping_arguments const & args, raptor_index<index_structure::hibf> & index);

} // namespace raptor
