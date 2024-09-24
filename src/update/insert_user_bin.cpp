// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/robin_hood.hpp>
#include <hibf/misc/empty_bins_by_fraction.hpp> // DEBUG

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>
#include <raptor/update/dump_index.hpp> // DEBUG

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
// static size_t depth = 0;
void get_ubs(robin_hood::unordered_flat_set<uint64_t> & ub_ids,
             raptor_index<index_structure::hibf> & index,
             size_t const ibf_idx)
{
    // ++depth;
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];
    // std::cerr << '[' << depth << "] user_bin_ids.size(): " << user_bin_ids.size() << '\n';
    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        size_t const ub = user_bin_ids[i];
        // std::cerr << '[' << depth << "] i: " << i << " UB: " << ub << '\n';
        // TODO: Should delete shrink?
        // -> Deleted bins don't have to be at the end!
        switch (ub)
        {
        case seqan::hibf::bin_kind::merged:
            // std::cerr << '[' << depth << "] Calling get_ubs for " << index.ibf().next_ibf_id[ibf_idx][i] << '\n';
            get_ubs(ub_ids, index, index.ibf().next_ibf_id[ibf_idx][i]);
            break;
        case seqan::hibf::bin_kind::deleted:
            break;
        default:
            // std::cerr << '[' << depth << "] Emplacing UB: " << ub << '\n';
            ub_ids.emplace(ub);
        }
    }
}

void get_overwrite_ibf(std::vector<size_t> & result, raptor_index<index_structure::hibf> & index, size_t const ibf_idx)
{
    result.push_back(ibf_idx);
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];
    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        if (user_bin_ids[i] == seqan::hibf::bin_kind::merged)
            get_overwrite_ibf(result, index, index.ibf().next_ibf_id[ibf_idx][i]);
    }
}

void partial_rebuild(update_arguments const & arguments,
                     detail::rebuild_location const & rebuild_location,
                     raptor_index<index_structure::hibf> & index)
{
    assert(index.ibf().ibf_bin_to_user_bin_id[rebuild_location.ibf_idx][rebuild_location.bin_idx]
           == seqan::hibf::bin_kind::merged);
    size_t const child_ibf_id = index.ibf().next_ibf_id[rebuild_location.ibf_idx][rebuild_location.bin_idx];

    std::vector<size_t> const ub_ids = [&]()
    {
        robin_hood::unordered_flat_set<uint64_t> ub_ids{};

        get_ubs(ub_ids, index, child_ibf_id);
        return std::vector<size_t>{ub_ids.begin(), ub_ids.end()};
    }();

    auto input_fn = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
    {
        raptor::file_reader<raptor::file_types::sequence> reader{index.shape(),
                                                                 static_cast<uint32_t>(index.window_size())};
        reader.hash_into(index.bin_path()[ub_ids[user_bin_id]], it);
    };
    // .number_of_hash_functions = ,
    // .relaxed_fpr = ,
    // .tmax =
    // TODO Store config in index
    seqan::hibf::config config{.input_fn = input_fn,
                               .number_of_user_bins = ub_ids.size(),
                               .maximum_fpr = index.fpr(),
                               .threads = arguments.threads,
                               .empty_bin_fraction = 0.85};
    auto config2 = config;
    config2.validate_and_set_defaults();
    std::cerr << "config2.tmax = " << config2.tmax << '\n';
    std::cerr << '\n';

    for (size_t i = 0; i < ub_ids.size(); ++i)
    {
        std::cerr << ub_ids[i] << ',';
    }
    std::cerr << '\n';

    seqan::hibf::hierarchical_interleaved_bloom_filter subindex{config};

    // DEV
    std::ranges::fill(subindex.ibf_vector[1].occupancy, 0u);
    std::ranges::fill(subindex.ibf_vector[1].occupied_bins, false);
    // ENDDEV
    std::cerr << "SUBINDEX\n";
    dump_index(subindex);

    // std::cout << subindex.ibf_vector.size() << '\n';
    // for (auto const & to_user_bin_id : subindex.ibf_bin_to_user_bin_id)
    // {
    //     std::cerr << "User bin[" << to_user_bin_id.size() << "]:   [";
    //     char sep{};
    //     for (auto const val : to_user_bin_id)
    //     {
    //         switch (val)
    //         {
    //         case seqan::hibf::bin_kind::deleted:
    //             std::cerr << sep << 'D';
    //             break;
    //         case seqan::hibf::bin_kind::merged:
    //             std::cerr << sep << 'M';
    //             break;
    //         default:
    //             std::cerr << sep << ub_ids[val];
    //         }

    //         sep = ',';
    //     }
    //     std::cerr << "]\n";
    // }
    // std::cerr << '\n';
    // std::cerr << "Number of IBFs: " << subindex.ibf_vector.size() << '\n';
    // std::cerr << "Next IBF ID:    ";
    // for (auto const val : subindex.next_ibf_id.front())
    //     std::cerr << val << ',';
    // std::cerr << "\nPrev IBF ID:    ";
    // for (auto const prev_pair : subindex.prev_ibf_id)
    //     std::cerr << prev_pair.ibf_idx << ':' << prev_pair.bin_idx << ',';
    // std::cerr << "\nIBF to UB:      ";
    // for (auto const val : subindex.ibf_bin_to_user_bin_id.front())
    //     std::cerr << val << ',';
    // std::cerr << "\n\n";

    auto & original_hibf = index.ibf();

    std::vector<size_t> overwrite_ibf_ids{};
    get_overwrite_ibf(overwrite_ibf_ids, index, child_ibf_id);
    std::cerr << "Overwrite IBF IDs: ";
    for (auto const val : overwrite_ibf_ids)
        std::cerr << val << ',';
    std::cerr << '\n';

    // TODO
    // While values in overwrite_ibf_ids, use these for new IBFs
    // If empty, use offset to push_back
    // If there are leftover values: ???
    // Or just delete and then append?
    // Maybe replace the first IBF, but delete its children and then append.
    // Children need to be modified anyway, but this way we don't need to modify the bookkeeping in the original too much.

    size_t const offset = original_hibf.ibf_vector.size() - 1u;
    // Handle the first IBF
    // Move ibf_vector
    original_hibf.ibf_vector[child_ibf_id] = std::move(subindex.ibf_vector[0]);
    // Move next_ibf_id
    auto & first_ibf_next_ibf_id = subindex.next_ibf_id[0];
    for (auto & id : first_ibf_next_ibf_id)
    {
        switch (id)
        {
        case 0:
            id = child_ibf_id;
            break;
        default:
            id += offset;
        }
    }
    original_hibf.next_ibf_id[child_ibf_id] = std::move(subindex.next_ibf_id[0]);
    // Move ibf_bin_to_user_bin_id
    auto & first_ibf_bin_to_user_bin_id = subindex.ibf_bin_to_user_bin_id[0];
    for (auto & id : first_ibf_bin_to_user_bin_id)
    {
        switch (id)
        {
        case seqan::hibf::bin_kind::deleted:
            break;
        case seqan::hibf::bin_kind::merged:
            break;
        default:
            id = ub_ids[id];
        }
    }
    original_hibf.ibf_bin_to_user_bin_id[child_ibf_id] = std::move(subindex.ibf_bin_to_user_bin_id[0]);
    // Prev_ibf_id does not change for first IBF
    // Move prev_ibf_id
    assert(subindex.ibf_vector[0].data() == nullptr);
    assert(subindex.next_ibf_id[0].empty());
    assert(subindex.ibf_bin_to_user_bin_id[0].empty());
    for (size_t i = 1; i < subindex.ibf_vector.size(); ++i)
    {
        original_hibf.ibf_vector.push_back(std::move(subindex.ibf_vector[i]));
        auto & ibf_next_ibf_id = subindex.next_ibf_id[i];
        for (auto & id : ibf_next_ibf_id)
        {
            id += offset;
        }
        original_hibf.next_ibf_id.push_back(std::move(subindex.next_ibf_id[i]));

        auto & ibf_bin_to_user_bin_id = subindex.ibf_bin_to_user_bin_id[i];
        for (auto & id : ibf_bin_to_user_bin_id)
        {
            switch (id)
            {
            case seqan::hibf::bin_kind::deleted:
                break;
            case seqan::hibf::bin_kind::merged:
                break;
            default:
                id = ub_ids[id];
            }
        }
        original_hibf.ibf_bin_to_user_bin_id.push_back(std::move(subindex.ibf_bin_to_user_bin_id[i]));

        auto prev_idx = subindex.prev_ibf_id[i];
        if (prev_idx.ibf_idx == 0)
            prev_idx.ibf_idx = child_ibf_id;
        else
            prev_idx.ibf_idx += offset;
        original_hibf.prev_ibf_id.push_back(prev_idx);
    }
}

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
    index.append_bin_path({arguments.user_bin_to_insert}); // TODO: update_bookkeeping, but it doesn't have the args
    auto const rebuild_location = detail::insert_tb_and_parents(kmers, insert_location, index);

    if (rebuild_location.ibf_idx != std::numeric_limits<size_t>::max())
    {
        std::cerr << "Partial Rebuild\n";
        partial_rebuild(arguments, rebuild_location, index);
        // Just construct an HIBF and attach it?
    }

    dump_index(index);

    //DEV
    partial_rebuild(arguments, detail::rebuild_location{0, 9}, index);

    // size_t const ibf_idx = std::get<0>(index_triple);
    // check if rebuilding is needed because of exceeding the tmax
    // check_tmax_rebuild(index, ibf_idx, update_arguments);
}

} // namespace raptor
