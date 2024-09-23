// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::insert_user_bin.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <hibf/contrib/robin_hood.hpp>

#include <raptor/argument_parsing/update_arguments.hpp>
#include <raptor/file_reader.hpp>
#include <raptor/index.hpp>

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

void get_ubs(robin_hood::unordered_flat_set<uint64_t> & ub_ids,
             raptor_index<index_structure::hibf> & index,
             size_t const ibf_idx)
{
    auto const & user_bin_ids = index.ibf().ibf_bin_to_user_bin_id[ibf_idx];
    for (size_t i = 0; i < user_bin_ids.size(); ++i)
    {
        size_t const ub = user_bin_ids[i];
        // TODO: Should delete shrink?
        // -> Deleted bins don't have to be at the end!
        switch (ub)
        {
        case seqan::hibf::bin_kind::merged:
            get_ubs(ub_ids, index, index.ibf().next_ibf_id[ibf_idx][i]);
            break;
        case seqan::hibf::bin_kind::deleted:
            break;
        default:
            ub_ids.emplace(ub);
        }
    }
}

void partial_rebuild(update_arguments const & arguments,
                     detail::rebuild_location const & rebuild_location,
                     raptor_index<index_structure::hibf> & index)
{
    assert(index.ibf().ibf_bin_to_user_bin_id[rebuild_location.ibf_idx][rebuild_location.bin_idx]
           == seqan::hibf::bin_kind::merged);
    size_t const child_ibf_id = index.ibf().next_ibf_id[rebuild_location.ibf_idx][rebuild_location.bin_idx];

    std::vector<size_t> ub_ids = [&]()
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
    seqan::hibf::config config{.input_fn = input_fn,
                               .number_of_user_bins = ub_ids.size(),
                               .maximum_fpr = index.fpr(),
                               .threads = arguments.threads,
                               .empty_bin_fraction = 0.0001};

    seqan::hibf::hierarchical_interleaved_bloom_filter subindex{config};

    std::cout << subindex.ibf_vector.size() << '\n';
    for (auto const & to_user_bin_id : subindex.ibf_bin_to_user_bin_id)
    {
        std::cerr << "User bin[" << to_user_bin_id.size() << "]:   [";
        char sep{};
        for (auto const val : to_user_bin_id)
        {
            switch (val)
            {
            case seqan::hibf::bin_kind::deleted:
                std::cerr << sep << 'D';
                break;
            case seqan::hibf::bin_kind::merged:
                std::cerr << sep << 'M';
                break;
            default:
                std::cerr << sep << ub_ids[val];
            }

            sep = ',';
        }
        std::cerr << "]\n";
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
    auto const rebuild_location = detail::insert_tb_and_parents(kmers, insert_location, index);

    if (rebuild_location.ibf_idx != std::numeric_limits<size_t>::max())
    {
        std::cerr << "Partial Rebuild\n";
        partial_rebuild(arguments, rebuild_location, index);
        // Just construct an HIBF and attach it?
    }

    //DEV
    partial_rebuild(arguments, detail::rebuild_location{0, 14}, index);

    // size_t const ibf_idx = std::get<0>(index_triple);
    // check if rebuilding is needed because of exceeding the tmax
    // check_tmax_rebuild(index, ibf_idx, update_arguments);
}

} // namespace raptor
