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
    dump_index(index.ibf());
}

void dump_index(seqan::hibf::hierarchical_interleaved_bloom_filter const & hibf)
{
    std::cerr << "\nDumping index\n";
    // std::cerr << "Window size: " << index.window_size() << '\n';
    // std::cerr << "Shape:       " << index.shape().to_string() << '\n';
    // std::cerr << "Parts:       " << index.parts() << '\n';
    // std::cerr << "FPR:         " << index.fpr() << '\n';
    // std::cerr << "Type:        " << (index.is_hibf() ? "HIBF" : "IBF") << '\n';
    // std::cerr << "Bin path:    " << index.bin_path().front().front() << '\n';
    for (size_t i = 0u; i < hibf.ibf_bin_to_user_bin_id.size(); ++i)
    {
        auto const & to_user_bin_id = hibf.ibf_bin_to_user_bin_id[i];

        std::cerr << "User bin " << i << " [" << to_user_bin_id.size() << "]:   [";
        char sep{};

        for (size_t j = 0u; j < to_user_bin_id.size(); ++j)
        {
            switch (to_user_bin_id[j])
            {
            case seqan::hibf::bin_kind::deleted:
                std::cerr << sep << 'D';
                break;
            case seqan::hibf::bin_kind::merged:
                std::cerr << sep << 'M' << hibf.next_ibf_id[i][j];
                break;
            default:
                std::cerr << sep << to_user_bin_id[j];
            }

            sep = ',';
        }
        std::cerr << "]\n";
    }

    size_t counter{};
    for (auto const & ibf : hibf.ibf_vector)
    {
        std::cerr << "IBF " << counter++ << '\n';
        std::cerr << "  Num bins: " << ibf.bin_count() << '\n';
        std::cerr << "  Occupancy[" << ibf.occupancy.size() << "]:   [";
        char sep{};
        for (auto const val : ibf.occupancy)
        {
            std::cerr << sep << val;
            sep = ',';
        }
        std::cerr << "]\n";
    }
    std::cerr << std::flush;
}

} // namespace raptor
