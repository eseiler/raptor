// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/misc/divide_and_ceil.hpp>

#include <raptor/index.hpp>

#include "strong_types.hpp"

namespace raptor::detail
{

// raptor::hierarchical_interleaved_bloom_filter::number_of_bins(size_t ibf_idx, size_t kmer_count)
size_t required_technical_bins(required_technical_bins_parameters const & params)
{
    auto compute_fpr = [&](size_t const elements)
    {
        double const exp_arg = (params.hash_count * elements) / static_cast<double>(params.bin_size);
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(params.hash_count * std::log(log_arg));
    };

    auto compute_split_fpr = [&](size_t const split)
    {
        double const fpr_tb = compute_fpr(seqan::hibf::divide_and_ceil(params.elements, split));
        return 1.0 - std::exp(std::log1p(-fpr_tb) * split);
    };

    size_t number_of_bins = seqan::hibf::divide_and_ceil(params.elements, params.max_elements);

    while (compute_split_fpr(number_of_bins) > params.fpr)
        ++number_of_bins;

    return number_of_bins;
}

size_t
find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t const ibf_idx, size_t const number_of_bins)
{
    // TODO: Increase more if empty_bin_fraction not satisfied
    // [[maybe_unused]] static constexpr double empty_bin_fraction = 0.0001; //TODO store in index
    auto & ibf = index.ibf().ibf_vector[ibf_idx];
    size_t const ibf_bin_count = [&]() -> size_t
    {
        auto search_result = std::ranges::search_n(ibf.occupancy, number_of_bins, 0u);
        if (!search_result.empty())
            return std::ranges::distance(ibf.occupancy.begin(), search_result.begin());
        return ibf.bin_count();
    }();

    // If nothing has been returned, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    size_t const new_bin_count{ibf_bin_count + number_of_bins};
    size_t const orig = ibf.bin_count();
    // If we can increase the number of bins without resizing the underlying bitvector
    if (ibf.try_increase_bin_number_to(seqan::hibf::bin_count{new_bin_count}))
    {
        // std::cerr << "[DEBUG] Try increase successful[" << ibf_idx << "]: " << orig << " to " << new_bin_count << '\n';
        // std::cerr << "[DEBUG] ibf_bin_count[" << ibf_idx << "]: " << ibf_bin_count << '\n';
        // std::cerr << "[DEBUG] number_of_bins[" << ibf_idx << "]: " << number_of_bins << '\n';
        // std::cerr << "[DEBUG] now[" << ibf_idx << "]: " << ibf.bin_count() << '\n';
        return ibf_bin_count;
    }
    {
        // std::cerr << "[DEBUG] Try increase NOT successful\n";
    }

    return std::numeric_limits<size_t>::max();
}

ibf_location find_ibf_size_splitting(std::vector<ibf_max> const & max_ibf_sizes,
                                     size_t const kmer_count,
                                     raptor_index<index_structure::hibf> & index)
{
    size_t const number_of_ibfs = max_ibf_sizes.size();

    // 1. Find the first IBF that is large enough to hold the new user bin.
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

        if (size_t const bin_idx = find_empty_bin_idx(index, ibf_idx, number_of_bins); bin_idx != std::numeric_limits<size_t>::max())
        {
            size_t const parent_ibf_idx = index.ibf().prev_ibf_id[ibf_idx].ibf_idx;
            // std::cerr << "[DEBUG] Case 1: " << ibf_idx << ' ' << parent_ibf_idx << '\n';
            return {.ibf_idx = ibf_idx, .bin_idx = bin_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
        }
    }
    while (ibf_size_idx-- != 0u);

    // 3. Check parent IBF and IBFs with sizes between first-large-enough IBF's size and parent's size.
    ibf_size_idx = binary_search_index;
    for (size_t ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx; ibf_idx != 0u;)
    {
        size_t const parent_ibf_idx = index.ibf().prev_ibf_id[ibf_idx].ibf_idx;
        auto const & ibf_parent = index.ibf().ibf_vector[parent_ibf_idx];

        // Parent has space.
        // this should trigger a resize and partial rebuild down the stream.
        if (parent_ibf_idx != 0)
        {
            if (size_t const bin_idx = find_empty_bin_idx(index, parent_ibf_idx, 1); bin_idx != std::numeric_limits<size_t>::max())
            {
                // this triggers resizing of the same ibf and leads to exceeding fpr -> full rebuild
                // return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
                // std::cerr << "[DEBUG] Case 2: " << ibf_idx << '\n';
                auto result = std::ranges::find_if(max_ibf_sizes, [parent_ibf_idx](ibf_max const & m) { return m.ibf_idx == parent_ibf_idx; });
                return {.ibf_idx = parent_ibf_idx, .bin_idx = bin_idx, .max_elements = max_ibf_sizes[std::ranges::distance(max_ibf_sizes.begin(), result)].max_elements};
            }
        }

        // Check IBFs with sizes between original IBF's size and parent's size.
        for (; ibf_size_idx < number_of_ibfs && max_ibf_sizes[ibf_size_idx].max_elements < ibf_parent.bin_size();
             ++ibf_size_idx)
        {
            ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
            if (size_t const bin_idx = find_empty_bin_idx(index, ibf_idx, 1); bin_idx != std::numeric_limits<size_t>::max())
            {
                // std::cerr << "[DEBUG] Case 3: " << ibf_idx << '\n';
                return {.ibf_idx = ibf_idx, .bin_idx = bin_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
            }
        }

        if (ibf_size_idx == number_of_ibfs)
            break;

        // if (ibf_size_idx == number_of_ibfs)
        // {
        //     // std::cerr << "[DEBUG] Case 4: " << ibf_idx << '\n';
        //     return {.ibf_idx = ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
        // }

        ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
    }

    // Any
    #if 1
    ibf_size_idx = binary_search_index;
    for (; ibf_size_idx < number_of_ibfs; ++ibf_size_idx)
    {
        size_t const ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx;
        if (size_t const bin_idx = find_empty_bin_idx(index, ibf_idx, 1); bin_idx != std::numeric_limits<size_t>::max())
        {
            // std::cerr << "[DEBUG] Case 5: " << ibf_idx << '\n';
            return {.ibf_idx = ibf_idx, .bin_idx = bin_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
        }
    }
    #endif

    // if (ibf_size_idx == number_of_ibfs)
    // {
    //     // std::cerr << "[DEBUG] Case 6: " << max_ibf_sizes[ibf_size_idx].ibf_idx << '\n';
    //     return {.ibf_idx = max_ibf_sizes[ibf_size_idx].ibf_idx, .max_elements = max_ibf_sizes[ibf_size_idx].max_elements};
    // }

    // a full rebuild will be triggered.
    // std::cerr << "[DEBUG] Case 7: 0\n";
    auto result = std::ranges::find_if(max_ibf_sizes, [](ibf_max const & m) { return m.ibf_idx == 0; });
    uint64_t bin_idx = find_empty_bin_idx(index, 0, 1);
    return {.ibf_idx = 0, .bin_idx = bin_idx, .max_elements = max_ibf_sizes[std::ranges::distance(max_ibf_sizes.begin(), result)].max_elements};
}

void update_bookkeeping(bookkeeping_arguments const & args, raptor_index<index_structure::hibf> & index)
{
    auto & ibf = index.ibf().ibf_vector[args.ibf_idx];
    size_t const new_number_of_bins = args.old_number_of_bins + args.number_of_new_bins;
    for (size_t i = args.old_number_of_bins; i < new_number_of_bins; ++i)
    {
        ibf.occupancy[i] = 1u;
    }
    index.ibf().next_ibf_id[args.ibf_idx].resize(new_number_of_bins, args.ibf_idx);
    index.ibf().ibf_bin_to_user_bin_id[args.ibf_idx].resize(new_number_of_bins, index.ibf().number_of_user_bins);
    index.ibf().number_of_user_bins += 1;
}

insert_location get_location(std::vector<ibf_max> const & max_ibf_sizes,
                             size_t const kmer_count,
                             raptor_index<index_structure::hibf> & index)
{
    auto const [ibf_idx, bidx, max_elements] = find_ibf_size_splitting(max_ibf_sizes, kmer_count, index);
    auto bin_idx = bidx;

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

    // uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, number_of_bins);

    // TODO: empty bins percentage
    // The current solution resizes the IBF here, but it would be more efficient to check the tmax function after calculating the new bin count and perhaps trigger a partial rebuild dirctly.
    if (bin_idx == std::numeric_limits<size_t>::max())
    {
        bin_idx = ibf.bin_count();
        // std::cerr << "[DEBUG] Forced increase[" << ibf_idx << "]: " << bin_idx << " to " << bin_idx + number_of_bins << '\n';
        ibf.increase_bin_number_to(seqan::hibf::bin_count{bin_idx + number_of_bins});
    }

    update_bookkeeping({.ibf_idx = ibf_idx, .old_number_of_bins = bin_idx, .number_of_new_bins = number_of_bins},
                       index);

    return insert_location{.ibf_idx = ibf_idx, .bin_idx = bin_idx, .number_of_bins = number_of_bins};
}

} // namespace raptor::detail
