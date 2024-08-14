// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <chrono>
#include <time.h>
#include "../../../include/cli_test.hpp"

struct timespec start, end;

struct update : public raptor_base, public testing::WithParamInterface<std::tuple<size_t, size_t, size_t>>
{};
// QUESTION how to measure RAM usage/peak memory
//1. Insert single bin and measure insertion time.
void insert_ub(std::string filename_ub, std::string filename_index){
    // 1.1 write filename to a file
    // in the location
    std::string insertion_paths = "data/tmp/insertion_paths";

    clock_gettime(CLOCK_MONOTONIC, &start); // 1.2 timer

    // TODO QUESTION how can i use execute app?
    // TEST_P(update, ebs_and_no_ebs){
    //cli_test_result const result_no_ebs =
//    cli_test_result execute_app
//    execute_app("raptor",
//                        "search",
//                        "--fpr 0.05",
//                        "--output result_no_ebs.out",
//                        "--hibf",
//                        "--index ",
//                        data(filename_index),
//                        "--bins ",
//                        insertion_paths;
//    EXPECT_EQ(result_no_ebs.out, std::string{});
//    EXPECT_EQ(result_no_ebs.err, std::string{});
//    RAPTOR_ASSERT_ZERO_EXIT(result_no_ebs);
    system("/usr/bin/time echo hello -f\"%M\" > cmd_test21.txt");

    clock_gettime(CLOCK_MONOTONIC, &end);

    // 1.3 calculate runtime
    double time_insertion = ((end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec)) * 1e-9;

}

void query_all_ubs(std::string filename_ub, std::string filename_index){
    // 2. Measure query time, averaged over all bins. Repeat.
    "series_of_insertions.fasta"
}

// 3. Save results to a file.

// I have not found a good way to track memory usage. We can just take the index size.  in CPP https://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-runtime-using-c, or call a command line tool outside.

//////////////////////////////////////////////////



TEST_P(update, ebs_and_no_ebs)
{

    auto const [number_of_repeated_bins, window_size, number_of_errors] = GetParam();

    // Without any empty bins
    cli_test_result const result_no_ebs =
        execute_app("raptor",
                    "search",
                    "--fpr 0.05",
                    "--output result_no_ebs.out",
                    "--error ",
                    std::to_string(number_of_errors),
                    "--p_max 0.4",
                    "--hibf",
                    "--index ",
                    data(hibf_no_ebs),
                    "--query ",
                    data(query_file));
    EXPECT_EQ(result_no_ebs.out, std::string{});
    EXPECT_EQ(result_no_ebs.err, std::string{});
    RAPTOR_ASSERT_ZERO_EXIT(result_no_ebs);
}


INSTANTIATE_TEST_SUITE_P(update_suite,
                         update,
                         testing::Combine(testing::Values(0, 16, 32), testing::Values(19, 23), testing::Values(0, 1)), // form parameters for different tests.
                         // --kmer 20, -- window 23 => parameters used for building
                         [](testing::TestParamInfo<update::ParamType> const & info)
                         {
                             std::string name = std::to_string(std::max<int>(1, std::get<0>(info.param) * 4)) + "_bins_"
                                              + std::to_string(std::get<1>(info.param)) + "_window_"
                                              + std::to_string(std::get<2>(info.param)) + "_error";
                             return name;
                         });