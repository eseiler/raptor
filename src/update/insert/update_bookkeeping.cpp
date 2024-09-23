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

// TODO make sure that values are actually written. occupancy is currently multiple of 64.
void update_bookkeeping(bookkeeping_arguments const & args, raptor_index<index_structure::hibf> & index)
{
    auto & ibf = index.ibf().ibf_vector[args.ibf_idx];
    // TODO: Should set_bin_number and increase_bin_number_to also handle occupancy (setting 0/true)?
    // increase_bin_number_to could, set_bin_number is also used for deletion (decrease number)
    size_t const new_number_of_bins = args.old_number_of_bins + args.number_of_new_bins;
    for (size_t i = args.old_number_of_bins; i < new_number_of_bins; ++i)
    {
        ibf.occupancy[i] = 0u;
        ibf.occupied_bins[i] = true;
    }
    index.ibf().next_ibf_id[args.ibf_idx].resize(new_number_of_bins, args.ibf_idx);
    index.ibf().ibf_bin_to_user_bin_id[args.ibf_idx].resize(new_number_of_bins, index.ibf().number_of_user_bins);
    index.ibf().number_of_user_bins += 1;
}

} // namespace raptor
