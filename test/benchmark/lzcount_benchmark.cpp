// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <random>

#include <raptor/test/benchmark/counting_vector.hpp>

static constexpr size_t operator""_MiB(unsigned long long int number)
{
    return number << 23;
}

static constexpr size_t vector_size{16384};

inline void fill_bitvector(raptor::benchmark::binning_bitvector & bitvector, size_t const bit_count)
{
    if (bit_count > 64)
        throw std::invalid_argument{"i must be in [0,64]"};

    static std::random_device rd{};
    static std::mt19937_64 gen{rd()};

    sdsl::bit_vector bits(64);
    for (size_t i = 0; i < bit_count; ++i)
        bits[i] = true;

    switch (bit_count)
    {
    case 0:
        std::fill(bitvector.begin(), bitvector.end(), 0);
        break;
    case 1:
        std::fill(bitvector.begin(), bitvector.end(), 1);
        break;
    default:
        for (size_t bit_pos = 0; bit_pos < vector_size; bit_pos += 64)
        {
            std::shuffle(bits.begin(), bits.end(), gen);
            bitvector.raw_data().set_int(bit_pos, bits.get_int(0));
        }
        break;
    }
}

template <typename counting_vector_t>
void run(benchmark::State & state)
{
    auto const bit_count = state.range(0);

    raptor::benchmark::binning_bitvector bitvector(vector_size);
    fill_bitvector(bitvector, bit_count);

    counting_vector_t counting_vector(vector_size);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(counting_vector += bitvector);
        benchmark::ClobberMemory();
    }
}

// tipping point: 36. >=37 invert_all, <=36 default
// 16384: 34

BENCHMARK_TEMPLATE(run, raptor::benchmark::counting_vector<uint16_t>)
    ->Arg(40)
    ->Arg(39)
    ->Arg(38)
    ->Arg(37)
    ->Arg(36)
    ->Arg(35)
    ->Arg(34)
    ->Arg(33)
    ->Arg(32);
// BENCHMARK_TEMPLATE(run, raptor::benchmark::counting_vector_invert<uint16_t>)
//     ->Arg(60)
//     ->Arg(50)
//     ->Arg(40)
//     ->Arg(30)
//     ->Arg(20)
//     ->Arg(10);
BENCHMARK_TEMPLATE(run, raptor::benchmark::counting_vector_invert_all<uint16_t>)
    ->Arg(40)
    ->Arg(39)
    ->Arg(38)
    ->Arg(37)
    ->Arg(36)
    ->Arg(35)
    ->Arg(34)
    ->Arg(33)
    ->Arg(32);

BENCHMARK_MAIN();
