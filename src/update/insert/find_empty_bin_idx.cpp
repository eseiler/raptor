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

size_t
find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t const ibf_idx, size_t const number_of_bins)
{
    // TODO: Increase more if empty_bin_fraction not satisfied
    [[maybe_unused]] static constexpr double empty_bin_fraction = 0.0001; //TODO store in index
    auto & ibf = index.ibf().ibf_vector[ibf_idx];
    size_t const ibf_bin_count = [&]() -> size_t
    {
        auto search_result = std::ranges::search_n(ibf.occupied_bins, number_of_bins, false);
        if (!search_result.empty())
            return std::ranges::distance(ibf.occupied_bins.begin(), search_result.begin());
        return ibf.bin_count();
    }();

    // If nothing has been returned, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    size_t const new_bin_count{ibf_bin_count + number_of_bins};
    // If we can increase the number of bins without resizing the underlying bitvector
    if (ibf.set_bin_count(seqan::hibf::bin_count{new_bin_count}))
    {
        return ibf_bin_count;
    }

    return std::numeric_limits<size_t>::max();
}

} // namespace raptor
