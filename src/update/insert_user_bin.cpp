// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/misc/divide_and_ceil.hpp>

#include <raptor/file_reader.hpp>
#include <raptor/update/insert_user_bin.hpp>

namespace raptor
{

robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::filesystem::path const & ub_file,
                                                       raptor_index<index_structure::hibf> const & index)
{
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    raptor::file_reader<raptor::file_types::sequence> reader{index.shape(), static_cast<uint32_t>(index.window_size())};
    reader.hash_into(ub_file, std::inserter(kmers, kmers.begin()));
    return kmers;
}

#if 0

std::tuple<uint64_t, uint64_t, uint16_t> get_location(size_t kmer_count,
                                                      robin_hood::unordered_flat_set<size_t> & kmers,
                                                      raptor_index<index_structure::hibf> & index,
                                                      update_arguments const & update_arguments)
{
    size_t root_idx = 0;
    size_t ibf_idx = -1;

    ibf_idx = find_ibf_size_splitting(kmers.size(), index, update_arguments);

    // calculate number of user bins needed.
    size_t number_of_bins = 1;
    if (index.ibf().ibf_max_kmers(ibf_idx) < (size_t)kmer_count)
    {
        number_of_bins = index.ibf().number_of_bins(ibf_idx, (int)kmer_count);
    }

    uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, number_of_bins);
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();

    // The current solution resizes the IBF here, but it would be more efficient to check the tmax function after calculating the new bin count and perhaps trigger a partial rebuild dirctly.
    if (bin_idx == ibf_bin_count)
    {
        size_t new_ibf_bin_count = new_bin_count(number_of_bins, update_arguments.empty_bin_percentage, ibf_bin_count);
        index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
    }

    return {ibf_idx, bin_idx, number_of_bins};
}

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto const kmers = compute_kmers(arguments.user_bin_to_insert, index);
    size_t const kmer_count = kmers.size();

    // This is the full code.
    // As a first step: Just look for a single empty bin

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

struct max_elements_parameters
{
    double fpr{};
    size_t hash_count{};
    size_t bin_size{};
};

// ceil(BITS / (-HASH / log(1 - exp(log(FPR) / HASH))))
inline constexpr size_t max_elements(max_elements_parameters const & params)
{
    assert(params.hash_count > 0);
    assert(params.fpr > 0.0);
    assert(params.fpr < 1.0);

    double const numerator{params.bin_size * std::log(1 - std::exp(std::log(params.fpr) / params.hash_count))};
    double const denominator{-static_cast<double>(params.hash_count)};

    double const result{std::ceil(numerator / denominator)};
    return result;
}

struct ibf_max
{
    size_t max_elements;
    size_t ibf_idx;

    constexpr auto operator<=>(ibf_max const & other) const = default;
};

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

struct required_technical_bins_parameters
{
    size_t bin_size{};
    size_t elements{};
    double fpr{};
    size_t hash_count{};
    size_t max_elements{};
};

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
    [[maybe_unused]] static constexpr double empty_bin_fraction = 0.0001; //TODO store in index
    auto & ibf = index.ibf().ibf_vector[ibf_idx];
    size_t const ibf_bin_count = ibf.bin_count();

    auto result = std::ranges::search_n(ibf.occupied_bins, number_of_bins, false);
    if (!result.empty())
        return std::ranges::distance(ibf.occupied_bins.begin(), result.begin());

    // If nothing has been returned, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    size_t const new_bin_count{ibf_bin_count + number_of_bins};
    // If we can increase the number of bins without resizing the underlying bitvector
    if (ibf.set_bin_count(seqan::hibf::bin_count{new_bin_count}))
    {
        // TODO make sure that values are actually written. occupancy is currently multiple of 64.
        ibf.occupancy.resize(new_bin_count);
        ibf.occupied_bins.resize(new_bin_count);
        index.ibf().next_ibf_id[ibf_idx].resize(new_bin_count, ibf_idx);
        index.ibf().ibf_bin_to_user_bin_id[ibf_idx].resize(new_bin_count, seqan::hibf::bin_kind::deleted);
        return ibf_bin_count;
    }

    return std::numeric_limits<size_t>::max();
    // Question resize above tmax?
}

size_t find_ibf_size_splitting(std::vector<ibf_max> const & array,
                               size_t const kmer_count,
                               raptor_index<index_structure::hibf> & index)
{
    // 1. FINDING THE IBF WITH BEST SIZE
    // https://godbolt.org/z/8dbznss37
    size_t const binary_search_index = [&]()
    {
        auto lower = std::ranges::lower_bound(array, ibf_max{.max_elements = kmer_count, .ibf_idx = 0});
        // There was no IBF with a size large enough to fit the new user bin.
        if (lower == array.end())
            return array.size() - 1u;
        else
            return static_cast<size_t>(std::ranges::distance(array.begin(), lower));
    }();

    // 2. EMPTY BINS IN SMALLER IBFS
    // index in the size_array of the 'perfect' IBF, which is just large enough to fit the new user bin.
    size_t smaller_ibf_idx = binary_search_index;
    // TODO: do-while loop
    while (true)
    {
        std::cerr << "smaller_ibf_idx = " << smaller_ibf_idx << '\n';
        size_t const ibf_idx = array[smaller_ibf_idx].ibf_idx;
        auto & ibf = index.ibf().ibf_vector[ibf_idx];

        // calculate among how many bins we should split
        size_t const number_of_bins = required_technical_bins({.bin_size = ibf.bin_size(),
                                                               .elements = kmer_count,
                                                               .fpr = index.fpr(),
                                                               .hash_count = ibf.hash_function_count(),
                                                               .max_elements = array[smaller_ibf_idx].max_elements});
        std::cerr << "Required bins = " << number_of_bins << '\n';

        // if we do find empty bins, just return this IBF.
        if (find_empty_bin_idx(index, ibf_idx, number_of_bins) != std::numeric_limits<size_t>::max())
            return ibf_idx;
        else if (smaller_ibf_idx == 0u)
            break;
        else
            smaller_ibf_idx -= 1;
    }

    // 3. EMPTY BINS IN LARGER IBFs
    size_t larger_ibf_idx = binary_search_index;
    size_t ibf_idx = array[larger_ibf_idx].ibf_idx;

    while (ibf_idx) // while the ibf is not the root
    {
        std::cerr << "larger_ibf_idx = " << larger_ibf_idx << '\n';
        size_t const parent_ibf_idx = index.ibf().prev_ibf_id[ibf_idx].ibf_idx;
        auto const & ibf_parent = index.ibf().ibf_vector[parent_ibf_idx];

        // the parent has at least one empty bin, which should accomodate a partial rebuild.
        if (find_empty_bin_idx(index, parent_ibf_idx, 1) != std::numeric_limits<size_t>::max())
        {
            return ibf_idx; // this should trigger a resize and partial rebuild down the stream.
        }
        else
        {
            size_t const size_parent = ibf_parent.bin_size();
            // and thereby larger_ibf_idx will be larger than the maximum value of the array
            while (larger_ibf_idx < array.size() && array[larger_ibf_idx].max_elements < size_parent)
            {
                ibf_idx = array[larger_ibf_idx].ibf_idx;
                if (find_empty_bin_idx(index, ibf_idx, 1) != std::numeric_limits<size_t>::max())
                    return ibf_idx; // if we do find empty bins, just return this IBF.
                else
                    larger_ibf_idx += 1;
            }
        }
        // this also includes the case where std::get<1>(array[larger_ibf_idx]) == parent_ibf_idx, which should be the case if the parents always have larger IBF bin sizes.
        if (larger_ibf_idx < array.size())
            ibf_idx = array[larger_ibf_idx].ibf_idx;
        else
            return ibf_idx;
    }
    std::cerr << "fallback = " << ibf_idx << '\n';
    return ibf_idx; // a full rebuild will be triggered.
}

void insert_user_bin(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    auto const kmers = compute_kmers(arguments.user_bin_to_insert, index);
    [[maybe_unused]] size_t const kmer_count = kmers.size();

    std::vector<ibf_max> const max_kmers = max_ibf_sizes(index);
    assert(std::ranges::is_sorted(max_kmers));
    std::cerr << "  Max kmers:   [";
    char sep{};
    for (auto && [max_kmer, ibf_idx] : max_kmers)
    {
        std::cerr << sep << max_kmer;
        sep = ',';
        (void)ibf_idx;
    }
    std::cerr << "]\n";

    // std::cerr << "Find 1 empty bin in top-level = " << find_empty_bin_idx(index, 0, 1) << '\n';
    // std::cerr << "Find 2 empty bins in top-level = " << find_empty_bin_idx(index, 0, 2) << '\n';
    std::cerr << "Find ibf for " << kmer_count << " kmers: " << find_ibf_size_splitting(max_kmers, kmer_count, index)
              << '\n';
    index.ibf().ibf_vector[0].occupied_bins.back() = true;
    index.ibf().ibf_vector[0].set_bin_count(seqan::hibf::bin_count{64u});
    std::cerr << "Find ibf for " << kmer_count << " kmers: " << find_ibf_size_splitting(max_kmers, kmer_count, index)
              << '\n';

    (void)arguments;
}

} // namespace raptor
