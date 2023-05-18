// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

/*!\file
 * \brief Implements raptor::chopper_layout.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#include <chopper/configuration.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <chopper/sketch/estimate_kmer_counts.hpp>
#include <chopper/sketch/execute.hpp>

#include <raptor/argument_parsing/init_shared_meta.hpp>
#include <raptor/layout/raptor_layout.hpp>
#include <chopper/layout/insert_empty_bins.hpp>

namespace raptor
{

void chopper_layout(sharg::parser & parser)
{
    chopper::configuration config;
    set_up_parser(parser, config);
    init_shared_meta(parser);
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    parser.parse();
    config.disable_sketch_output = !parser.is_option_set("output-sketches-to");

    chopper::layout::layout hibf_layout{};
    std::vector<std::string> filenames{};
    std::vector<size_t> kmer_counts{};
    std::vector<chopper::sketch::hyperloglog> sketches{};
    std::vector<bool> empty_bins; // A bitvector indicating whether a bin is empty (1) or not (0).
    std::vector<size_t> empty_bin_cum_sizes; // The cumulative k-mer count for the first empty bin until empty bin i
    chopper::layout::insert_empty_bins(empty_bins, empty_bin_cum_sizes,
                      kmer_counts, sketches, filenames, config);

    chopper::sketch::read_data_file(config, filenames);

    chopper::sketch::execute(config, filenames, sketches);
    chopper::sketch::estimate_kmer_counts(sketches, kmer_counts);

    chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                              .hibf_layout = &hibf_layout,
                              .kmer_counts = kmer_counts,
                              .sketches = sketches,
                              .empty_bins = empty_bins,
                              .empty_bin_cum_sizes = empty_bin_cum_sizes};

    chopper::layout::execute(config, filenames, store);
}

}