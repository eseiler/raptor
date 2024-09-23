// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include "insert/fwd.hpp"

namespace raptor
{

#if 0

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto const kmers = compute_kmers(arguments.user_bin_to_insert, index);
    size_t const kmer_count = kmers.size();

    std::tuple<uint64_t, uint64_t, uint16_t> index_triple = get_location(kmer_count, kmers, index, update_arguments);
    std::tuple<uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers, index_triple, index);

    index.ibf().user_bins.update_filename_indices(filename_cast, index_triple);
    update_sketch(filename, update_arguments, index, false); //Compute and store sketches for the new UBs

    // If the ibf_idx is lower than the total number of IBFs (the initilization value), then the respective IBF needs to be rebuild.
    if (std::get<0>(rebuild_index_tuple) < index.ibf().ibf_vector.size())
        partial_rebuild(rebuild_index_tuple, index, update_arguments);

    size_t const ibf_idx = std::get<0>(index_triple);

    // check if rebuilding is needed because of exceeding the tmax

    check_tmax_rebuild(index, ibf_idx, update_arguments);
}

#endif

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto const kmers = compute_kmers(arguments.user_bin_to_insert, index);
    size_t const kmer_count = kmers.size();

    std::vector<ibf_max> const max_kmers = max_ibf_sizes(index);
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

    auto const insert_location = get_location(max_kmers, kmer_count, index);
    auto const rebuild_location = insert_tb_and_parents(kmers, insert_location, index);

    if (rebuild_location.ibf_idx != std::numeric_limits<size_t>::max())
    {
        std::cerr << "Partial Rebuild\n";
        // partial_rebuild(rebuild_location, index);
    }
}

} // namespace raptor
