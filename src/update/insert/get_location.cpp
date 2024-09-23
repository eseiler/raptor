// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include "fwd.hpp"

namespace raptor
{

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index)
{
    auto const [ibf_idx, max_elements] = find_ibf_size_splitting(max_ibf_sizes, kmer_count, index);

    auto & ibf = index.ibf().ibf_vector[ibf_idx];

    // calculate number of user bins needed.
    size_t number_of_bins = 1;
    if (max_elements < kmer_count)
    {
        number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                  .elements = kmer_count,
                                                  .fpr = index.fpr(),
                                                  .hash_count = ibf.hash_function_count(),
                                                  .max_elements = max_elements});
    }

    uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, number_of_bins);

    // TODO: empty bins percentage
    // The current solution resizes the IBF here, but it would be more efficient to check the tmax function after calculating the new bin count and perhaps trigger a partial rebuild dirctly.
    if (bin_idx == std::numeric_limits<size_t>::max())
    {
        bin_idx = ibf.bin_count();
        ibf.increase_bin_number_to(seqan::hibf::bin_count{bin_idx + number_of_bins});
    }

    update_bookkeeping({.ibf_idx = ibf_idx, .old_number_of_bins = bin_idx, .number_of_new_bins = number_of_bins},
                       index);

    return insert_location{.ibf_idx = ibf_idx, .bin_idx = bin_idx, .number_of_bins = number_of_bins};
}

} // namespace raptor
