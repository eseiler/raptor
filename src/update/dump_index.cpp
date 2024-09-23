// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Implements raptor::dump_index.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <raptor/update/dump_index.hpp>

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
                std::cerr << sep << val;
            }

            sep = ',';
        }
        std::cerr << "]\n";
    }
    for (auto const & ibf : index.ibf().ibf_vector)
    {
        std::cerr << "Num bins:    " << ibf.bin_count() << '\n';
        std::cerr << "  Occupancy[" << ibf.occupancy.size() << "]:   [";
        char sep{};
        for (auto const val : ibf.occupancy)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
        std::cerr << "  Occupied[" << ibf.occupied_bins.size() << "]:    [";
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

} // namespace raptor
