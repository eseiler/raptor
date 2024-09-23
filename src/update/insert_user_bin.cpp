// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

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

robin_hood::unordered_flat_set<uint64_t> compute_kmers(std::filesystem::path const & ub_file,
                                                       raptor_index<index_structure::hibf> const & index)
{
    robin_hood::unordered_flat_set<uint64_t> kmers{};
    raptor::file_reader<raptor::file_types::sequence> reader{index.shape(), static_cast<uint32_t>(index.window_size())};
    reader.hash_into(ub_file, std::inserter(kmers, kmers.begin()));
    return kmers;
}

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

struct ibf_location
{
    size_t ibf_idx;
    size_t max_elements;
};

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

    for (size_t i = insert_location.bin_idx; i < insert_location.bin_idx + insert_location.number_of_bins; ++i)
    {
        if (compute_fpr(ibf, i) > index.fpr())
        {
            std::cerr << "Rebuild needed" << '\n';
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
    rebuild_location rebuild_location{.ibf_idx = index.ibf().ibf_vector.size(), .bin_idx = 0};
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
    [[maybe_unused]] auto const rebuild_location = insert_tb_and_parents(kmers, insert_location, index);

    (void)arguments;
}

} // namespace raptor
