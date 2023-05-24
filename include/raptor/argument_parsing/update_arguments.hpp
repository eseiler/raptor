// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <filesystem>
#include <vector>

#include <seqan3/search/kmer_index/shape.hpp>

namespace raptor
{
/*!\brief Arguments required for the Dynamic HIBF.
* \author Myrthe Willemsen
*/
struct update_arguments
{
    // needed for k-mer loading
    seqan3::shape shape{seqan3::ungapped{20}};
    uint32_t window_size{23};



    uint8_t parts{1u};
    //!\brief The index is compressed
    bool compressed{false};
    //!\brief Bin files are extracted from the `bin_file` and stored here.
    std::vector<std::vector<std::string>> bin_path{};
    bool is_minimiser{false};
    // Note, the above: `window_size`, `compressed`, `parts`,  `bin_path` and `shape` stored as part of the `index` datastructure. Adding them as input again is redundant. `kmer_size` could be obtained from `shape`,

    //!\brief The index is an HIBF
    bool is_hibf{false};

    //!\brief The method to select the IBF when inserting new user bins: "find_ibf_idx_traverse_by_similarity"; "find_ibf_idx_ibf_size"; "find_ibf_idx_traverse_by_fpr";
    std::string ibf_selection_method{"find_ibf_idx_traverse_by_fpr"};

    //!\brief The percentage of empty bins sampled during layout computation.
    double empty_bin_percentage{0.1};

    //!\brief Should updates account for sequence similarities, similar to setting the `rearrange_user_bins` argument in the layout algorithm.
    bool similarity{false};

    //!\brief Maximum number of technical bins per IBF. Needed for layouting.
    double tmax{64};
    bool tmax_condition{true}; //Indicates if the tmax should also be a restriction on IBFs at lower levels

    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits{12};

    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t threads{1u};


    // Following arguments are mutually exclusive and indicate what update operation should be performed
    bool delete_ubs{false};
    bool insert_ubs{false};
    bool insert_sequences{false};
    bool delete_sequences{false};

    // Filenames
    std::filesystem::path bin_file{};
    std::filesystem::path in_file{};
    std::filesystem::path out_file{};
    std::filesystem::path sketch_directory{"sketches"};
    std::string insert_sequence_appendix{"_insertsequences"};
};

} // namespace raptor
