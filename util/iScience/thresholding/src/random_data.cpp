// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Random data generator.
 */

#include <filesystem>
#include <random>

#include <sharg/all.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

void run_program(std::filesystem::path const & out_directory,
                 size_t const reference_length,
                 size_t const number_of_queries,
                 size_t const query_length,
                 uint8_t const min_error,
                 uint8_t const max_error,
                 size_t const seed)
{
    std::mt19937_64 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_dna(0, 3);                                  // dna4
    std::uniform_int_distribution<size_t> dis_pos(0, reference_length - query_length - 1); // positions
    std::uniform_int_distribution<size_t> dis_err(0, query_length - 1);                    // errors

    // =========================================================================
    // Reference
    // =========================================================================
    std::vector<seqan3::dna4> reference_sequence;
    reference_sequence.reserve(reference_length);

    for (size_t i = 0; i < reference_length; ++i)
        reference_sequence.push_back(seqan3::dna4{}.assign_rank(dis_dna(gen)));

    {
        std::filesystem::path reference_filename = out_directory / "reference.fasta";
        seqan3::sequence_file_output reference_out{reference_filename};
        reference_out.emplace_back(reference_sequence, "reference_" + std::to_string(reference_length));
    }

    // =========================================================================
    // Queries with fixed number of errors
    // =========================================================================
    std::vector<seqan3::dna4> query_sequence;
    query_sequence.reserve(query_length);

    for (uint8_t errors = min_error; errors <= max_error; ++errors)
    {
        std::filesystem::path query_filename = out_directory / ("query_e" + std::to_string(errors) + ".fasta");
        {
            seqan3::sequence_file_output query_out{query_filename};
            for (size_t i = 0; i < number_of_queries; ++i)
            {
                size_t start = dis_pos(gen);
                auto query_span = std::span{reference_sequence.data() + start, query_length};
                query_sequence.assign(query_span.begin(), query_span.end());

                for (size_t j = 0; j < errors; ++j)
                {
                    size_t error_pos = dis_err(gen);
                    seqan3::dna4 new_base{};
                    new_base.assign_rank(dis_dna(gen));
                    while (new_base == query_sequence[error_pos])
                        new_base.assign_rank(dis_dna(gen));
                    query_sequence[error_pos] = new_base;
                }
                query_out.emplace_back(query_sequence, "query_" + std::to_string(i));
            }
        }
    }

    // =========================================================================
    // Random queries (to check for false positives)
    // =========================================================================
    {
        std::filesystem::path query_filename = out_directory / "query_random.fasta";
        seqan3::sequence_file_output query_out{query_filename};

        for (size_t i = 0; i < number_of_queries; ++i)
        {
            query_sequence.clear();

            for (size_t i = 0; i < query_length; ++i)
                query_sequence.push_back(seqan3::dna4{}.assign_rank(dis_dna(gen)));

            query_out.emplace_back(query_sequence, "query_" + std::to_string(i));
        }
    }
}

struct cmd_arguments
{
    std::filesystem::path output_file{};
    size_t reference_length{10'000'000};
    size_t number_of_queries{100'000};
    size_t query_length{150};
    uint8_t min_error{0};
    uint8_t max_error{3};
    size_t seed{0x7E82B6F280D4706B};
};

void initialize_argument_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.short_description = "This programs creates random data for testing purposes.";
    parser.info.version = "1.0.0";

    parser.add_option(args.output_file,
                      sharg::config{.short_id = '\0',
                                    .long_id = "out",
                                    .description = "The output directory where the files will be located.",
                                    .required = true,
                                    .validator = sharg::output_directory_validator{}});
    parser.add_option(args.reference_length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "reference-size",
                                    .description = "The length of the reference.",
                                    .validator = sharg::arithmetic_range_validator{1, 1'000'000'000}});
    parser.add_option(args.number_of_queries,
                      sharg::config{.short_id = '\0',
                                    .long_id = "number-of-queries",
                                    .description = "The number of queries.",
                                    .validator = sharg::arithmetic_range_validator{1, 1'000'000'000}});
    parser.add_option(args.query_length,
                      sharg::config{.short_id = '\0',
                                    .long_id = "query_length",
                                    .description = "The length of the queries.",
                                    .validator = sharg::arithmetic_range_validator{1, 1'000'000}});
    parser.add_option(
        args.min_error,
        sharg::config{.short_id = '\0', .long_id = "min_error", .description = "The minimal number of errors."});
    parser.add_option(
        args.max_error,
        sharg::config{.short_id = '\0', .long_id = "max_error", .description = "The maximal number of errors."});
    parser.add_option(
        args.seed,
        sharg::config{.short_id = '\0', .long_id = "seed", .description = "The seed to use.", .advanced = true});
}

int main(int argc, char ** argv)
{
    sharg::parser myparser{"Minimizers", argc, argv, sharg::update_notifications::off};
    cmd_arguments args{};

    initialize_argument_parser(myparser, args);

    try
    {
        myparser.parse();
        if (args.min_error > args.max_error)
            throw sharg::user_input_error{"The minimum number of errors cannot be greater than the maximum."};
    }
    catch (sharg::parser_error const & ext)
    {
        std::cerr << "[Error] " << ext.what() << "\n";
        return -1;
    }

    run_program(args.output_file,
                args.reference_length,
                args.number_of_queries,
                args.query_length,
                args.min_error,
                args.max_error,
                args.seed);

    return 0;
}
