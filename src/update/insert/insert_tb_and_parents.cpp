// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/build/insert_into_ibf.hpp>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

void insert_into_ibf(robin_hood::unordered_flat_set<size_t> const & kmers,
                     insert_location const & insert_location,
                     raptor_index<index_structure::hibf> & index,
                     rebuild_location & rebuild_index_tuple)
{
    auto & ibf = index.ibf().ibf_vector[insert_location.ibf_idx];

    seqan::hibf::build::insert_into_ibf(kmers,
                                        insert_location.number_of_bins,
                                        insert_location.bin_idx,
                                        ibf,
                                        index.ibf().fill_ibf_timer);

    auto compute_fpr = [](auto const & ibf, size_t const bin_idx)
    {
        double const exp_arg =
            (ibf.hash_function_count() * ibf.occupancy[bin_idx]) / static_cast<double>(ibf.bin_size());
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(ibf.hash_function_count() * std::log(log_arg));
    };

    // TODO: Won't the kmers be evenly split? In this case, one computation is enough.
    for (size_t i = insert_location.bin_idx; i < insert_location.bin_idx + insert_location.number_of_bins; ++i)
    {
        if (compute_fpr(ibf, i) > index.fpr())
        {
            rebuild_index_tuple.ibf_idx = insert_location.ibf_idx;
            rebuild_index_tuple.bin_idx = insert_location.bin_idx;
            break;
        }
    }

    // TODO
    // // to improve the implementation, Perhaps do the FPR calculations for all bins to which kmers will be inserted before actually inserting.
    // index.ibf().update_occupancy_table(kmers.size()-union_count, ibf_idx, start_bin_idx, number_of_bins);
    // auto fpr = index.ibf().update_fpr(ibf_idx, start_bin_idx, number_of_bins); // this should be done after updating the occupancy table.

    // // Keep: If any fpr is too high, the index needs to be rebuild.
    // if (fpr > index.ibf().fpr_max){
    //     //assert(index.ibf().next_ibf_id[ibf_idx][start_bin_idx] != ibf_idx); //assert that it is not a leaf bin. fpr should not reach fpr max for the leaf bin ibf, because we search a location such that it is 'feasible', i.e. binsize should be sufficient to accomodate new UB. However, this is not the case if one also does sequence insertions.
    //     rebuild_index_tuple = std::make_tuple(ibf_idx, start_bin_idx);
    // }
}

rebuild_location insert_tb_and_parents(robin_hood::unordered_flat_set<size_t> const & kmers,
                                       insert_location insert_location,
                                       raptor_index<index_structure::hibf> & index)
{
    rebuild_location rebuild_location{.ibf_idx = std::numeric_limits<size_t>::max(),
                                      .bin_idx = std::numeric_limits<size_t>::max()};
    while (true)
    {
        insert_into_ibf(kmers, insert_location, index, rebuild_location);
        if (insert_location.ibf_idx == 0u)
            break;
        auto const parent = index.ibf().prev_ibf_id[insert_location.ibf_idx];
        insert_location.ibf_idx = parent.ibf_idx;
        insert_location.bin_idx = parent.bin_idx;
        insert_location.number_of_bins = 1u; // number of bins will be 1 for the merged bins (i assume)
    }

    return rebuild_location;
}

} // namespace raptor::detail
