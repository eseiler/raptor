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

ibf_location find_ibf_size_splitting(std::vector<ibf_max> const & max_ibf_sizes,
                                     size_t const kmer_count,
                                     raptor_index<index_structure::hibf> & index)
{
    size_t const number_of_ibfs = max_ibf_sizes.size();

    // 1. Find best fit.
    // https://godbolt.org/z/8dbznss37
    size_t const binary_search_index = [&]()
    {
        auto lower = std::ranges::lower_bound(max_ibf_sizes, ibf_max{.max_elements = kmer_count, .ibf_idx = 0});
        // There was no IBF with a size large enough to fit the new user bin.
        if (lower == max_ibf_sizes.end())
            return number_of_ibfs - 1u;
        else
            return static_cast<size_t>(std::ranges::distance(max_ibf_sizes.begin(), lower));
    }();

    // 2. Check smaller IBFs.
    size_t ibf_size_idx = binary_search_index;

    do
    {
        size_t const ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
        auto & ibf = index.ibf().ibf_vector[ibf_idx];

        size_t const number_of_bins =
            required_technical_bins({.bin_size = ibf.bin_size(),
                                     .elements = kmer_count,
                                     .fpr = index.fpr(),
                                     .hash_count = ibf.hash_function_count(),
                                     .max_elements = max_ibf_sizes[ibf_size_idx].max_elements});

        if (find_empty_bin_idx(index, ibf_idx, number_of_bins) != std::numeric_limits<size_t>::max())
            return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
    }
    while (ibf_size_idx-- != 0u);

    // 3. Check parent IBF and IBFs IBFs with sizes between original IBF's size and parent's size.
    ibf_size_idx = binary_search_index;
    for (size_t ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx; ibf_idx != 0u;)
    {
        size_t const parent_ibf_idx = index.ibf().prev_ibf_id[ibf_idx].ibf_idx;
        auto const & ibf_parent = index.ibf().ibf_vector[parent_ibf_idx];

        // Parent has space.
        // this should trigger a resize and partial rebuild down the stream.
        if (find_empty_bin_idx(index, parent_ibf_idx, 1) != std::numeric_limits<size_t>::max())
            return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};

        // Check IBFs with sizes between original IBF's size and parent's size.
        for (; ibf_size_idx < number_of_ibfs && max_ibf_sizes[ibf_size_idx].max_elements < ibf_parent.bin_size();
             ++ibf_size_idx)
        {
            ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
            if (find_empty_bin_idx(index, ibf_idx, 1) != std::numeric_limits<size_t>::max())
                return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
        }

        if (ibf_size_idx == number_of_ibfs)
            return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};

        ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
    }

    // a full rebuild will be triggered.
    return {.ibf_idx = 0, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
}

} // namespace raptor
