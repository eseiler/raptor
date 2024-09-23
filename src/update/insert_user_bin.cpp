// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/robin_hood.hpp>

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>

#include "insert/strong_types.hpp"

namespace raptor::detail
{

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index);

rebuild_location insert_tb_and_parents(robin_hood::unordered_flat_set<size_t> const & kmers,
                                       insert_location insert_location,
                                       raptor_index<index_structure::hibf> & index);

robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::filesystem::path const & ub_file,
                                                       raptor_index<index_structure::hibf> const & index)
{
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    raptor::file_reader<raptor::file_types::sequence> reader{index.shape(), static_cast<uint32_t>(index.window_size())};
    reader.hash_into(ub_file, std::inserter(kmers, kmers.begin()));
    return kmers;
}

// ceil(BITS / (-HASH / log(1 - exp(log(FPR) / HASH))))
size_t max_elements(max_elements_parameters const & params)
{
    assert(params.hash_count > 0);
    assert(params.fpr > 0.0);
    assert(params.fpr < 1.0);

    double const numerator{params.bin_size * std::log(1 - std::exp(std::log(params.fpr) / params.hash_count))};
    double const denominator{-static_cast<double>(params.hash_count)};

    double const result{std::ceil(numerator / denominator)};
    return result;
}

// Shouldn't this be the occupancy? However, bins might get cleared.
// Answer: Split bin correction :)
std::vector<ibf_max> max_ibf_sizes(raptor_index<index_structure::hibf> const & index)
{
    auto const & ibf_vector = index.ibf().ibf_vector;
    std::vector<ibf_max> max_sizes{};
    max_sizes.reserve(ibf_vector.size());

    for (size_t i = 0; i < ibf_vector.size(); ++i)
    {
        auto const & ibf = ibf_vector[i];
        size_t const max_kmers = max_elements({.fpr = index.fpr(), //
                                               .hash_count = ibf.hash_function_count(),
                                               .bin_size = ibf.bin_size()});
        max_sizes.push_back({.max_elements = max_kmers, .ibf_idx = i});
    }
    std::ranges::sort(max_sizes);
    return max_sizes;
}

} // namespace raptor::detail

namespace raptor
{

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto const kmers = detail::compute_kmers(arguments.user_bin_to_insert, index);
    size_t const kmer_count = kmers.size();

    std::vector<detail::ibf_max> const max_kmers = detail::max_ibf_sizes(index);
    assert(std::ranges::is_sorted(max_kmers));
    // std::cerr << "  Max kmers:   [";
    // char sep{};
    // for (auto && [max_kmer, ibf_idx] : max_kmers)
    // {
    //     std::cerr << sep << max_kmer;
    //     sep = ',';
    //     (void)ibf_idx;
    // }
    // std::cerr << "]\n";

    auto const insert_location = detail::get_location(max_kmers, kmer_count, index);
    auto const rebuild_location = detail::insert_tb_and_parents(kmers, insert_location, index);

    if (rebuild_location.ibf_idx != std::numeric_limits<size_t>::max())
    {
        std::cerr << "Partial Rebuild\n";
        // partial_rebuild(rebuild_location, index);
        // Just construct an HIBF and attach it?
    }

    // size_t const ibf_idx = std::get<0>(index_triple);
    // check if rebuilding is needed because of exceeding the tmax
    // check_tmax_rebuild(index, ibf_idx, update_arguments);
}

} // namespace raptor
