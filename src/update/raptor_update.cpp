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

template <typename index_t>
void dump_index(std::filesystem::path const & index_filename) = delete;

template <>
void dump_index<index_structure::ibf>(std::filesystem::path const & index_filename)
{
    std::ifstream index_file{index_filename};
    cereal::BinaryInputArchive archive{index_file};
    raptor::raptor_index<index_structure::ibf> index;
    archive(index);
    auto & ibf = index.ibf();
    std::cerr << ibf.bin_count() << '\n';
    for (auto val : ibf.occupancy)
    {
        std::cout << val << ',';
    }
    std::cout << '\n';
    for (auto val : ibf.occupied_bins)
    {
        std::cout << val << ',';
    }
    std::cout << '\n';
}

template <>
void dump_index<index_structure::hibf>(std::filesystem::path const & index_filename)
{
    std::ifstream index_file{index_filename};
    cereal::BinaryInputArchive archive{index_file};
    raptor::raptor_index<index_structure::hibf> index;
    archive(index);
    for (auto & ibf : index.ibf().ibf_vector)
    {
        std::cerr << ibf.bin_count() << '\n';
        for (auto val : ibf.occupancy)
        {
            std::cout << val << ',';
        }
        std::cout << '\n';
        for (auto val : ibf.occupied_bins)
        {
            std::cout << val << ',';
        }
        std::cout << '\n';
    }
}

void raptor_update(update_arguments const & arguments)
{
    if (arguments.is_hibf)
        dump_index<index_structure::hibf>(arguments.index_file);
    else
        dump_index<index_structure::ibf>(arguments.index_file);

    return;
}

} // namespace raptor
