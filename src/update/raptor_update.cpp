// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::raptor_update.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <fstream>

#include <raptor/index.hpp>
#include <raptor/update/update.hpp>

namespace raptor
{

void dump_index(raptor_index<index_structure::hibf> const & index)
{
    std::cerr << "\nDumping index\n";
    // std::cerr << "Window size: " << index.window_size() << '\n';
    // std::cerr << "Shape:       " << index.shape().to_string() << '\n';
    // std::cerr << "Parts:       " << index.parts() << '\n';
    // std::cerr << "FPR:         " << index.fpr() << '\n';
    // std::cerr << "Type:        " << (index.is_hibf() ? "HIBF" : "IBF") << '\n';
    // std::cerr << "Bin path:    " << index.bin_path().front().front() << '\n';
    for (auto const & to_user_bin_id : index.ibf().ibf_bin_to_user_bin_id)
    {
        std::cerr << "User bin:   [";
        char sep{};
        for (auto const val : to_user_bin_id)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
    }
    for (auto const & ibf : index.ibf().ibf_vector)
    {
        std::cerr << "Num bins:    " << ibf.bin_count() << '\n';
        std::cerr << "  Occupancy:   [";
        char sep{};
        for (auto const val : ibf.occupancy)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
        std::cerr << "  Occupied:    [";
        sep = {};
        for (auto const val : ibf.occupied_bins)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
    }
    std::cerr << std::flush;
}

void delete_user_bins(update_arguments const & arguments, raptor_index<index_structure::hibf> & index)
{
    std::cerr << "\nDeleting user bins: ";
    for (auto const bin : arguments.user_bins_to_delete)
        std::cerr << bin << ' ';
    std::cerr << '\n';

    auto & ibf_bin_to_user_bin_id = index.ibf().ibf_bin_to_user_bin_id;
    auto & next_ibf_id = index.ibf().next_ibf_id;
    std::vector<seqan::hibf::bin_index> technical_bins_to_delete{};

    for (size_t ibf_index = 0; ibf_index < ibf_bin_to_user_bin_id.size(); ++ibf_index)
    {
        technical_bins_to_delete.clear();

        auto & to_user_bin_id = ibf_bin_to_user_bin_id[ibf_index];
        auto & ibf = index.ibf().ibf_vector[ibf_index];

        for (size_t bin_index = 0; bin_index < to_user_bin_id.size(); ++bin_index)
        {
            int64_t & user_bin_id = to_user_bin_id[bin_index];
            if (std::ranges::find(arguments.user_bins_to_delete, user_bin_id) != arguments.user_bins_to_delete.end())
            {
                technical_bins_to_delete.push_back(seqan::hibf::bin_index{bin_index});
                user_bin_id = -2;
            }
        }

        if (!technical_bins_to_delete.empty())
        {
            ibf.clear(technical_bins_to_delete);
            for (auto const technical_bin_index : technical_bins_to_delete)
            {
                ibf.occupancy[technical_bin_index.value] = 0u;
                ibf.occupied_bins[technical_bin_index.value] = false;
            }
            ibf.overwrite_bin_count(ibf.bin_count() - technical_bins_to_delete.size());
            assert(ibf.occupied_bins.none() ^ ibf.bin_count());
            // Delete in parent
            if (ibf_index != 0 && ibf.occupied_bins.none())
            {
                for (size_t parent_ibf_idx = 0; parent_ibf_idx < next_ibf_id.size(); ++parent_ibf_idx)
                {
                    if (parent_ibf_idx == ibf_index)
                        continue;

                    auto & parent_next_ibf = next_ibf_id[parent_ibf_idx];

                    for (size_t parent_bin_index = 0; parent_bin_index < parent_next_ibf.size(); ++parent_bin_index)
                    {
                        if (parent_next_ibf[parent_bin_index] == static_cast<int64_t>(ibf_index))
                        {
                            auto & parent_ibf = index.ibf().ibf_vector[parent_ibf_idx];
                            parent_ibf.clear(seqan::hibf::bin_index{parent_bin_index});

                            parent_ibf.occupancy[parent_bin_index] = 0u;
                            parent_ibf.occupied_bins[parent_bin_index] = false;
                            ibf_bin_to_user_bin_id[parent_ibf_idx][parent_bin_index] = -2;
                            parent_ibf.overwrite_bin_count(parent_ibf.bin_count() - 1u);
                        }
                    }
                }
            }
        }
    }
}

void raptor_update(update_arguments const & arguments)
{
    std::ifstream index_file{arguments.index_file};
    cereal::BinaryInputArchive archive{index_file};
    raptor::raptor_index<index_structure::hibf> index;
    archive(index);

    dump_index(index);
    if (!arguments.user_bins_to_delete.empty())
    {
        delete_user_bins(arguments, index);
        dump_index(index);
    }
}

} // namespace raptor
