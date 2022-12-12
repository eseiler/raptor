// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <sdsl/int_vector.hpp>

namespace raptor::benchmark
{

class binning_bitvector
{
private:
    using data_type = sdsl::bit_vector;
    data_type data{};

public:
    binning_bitvector() = default;
    binning_bitvector(binning_bitvector const &) = default;
    binning_bitvector & operator=(binning_bitvector const &) = default;
    binning_bitvector(binning_bitvector &&) = default;
    binning_bitvector & operator=(binning_bitvector &&) = default;
    ~binning_bitvector() = default;

    explicit binning_bitvector(size_t const size) : data(size)
    {}

    size_t size() const noexcept
    {
        return data.size();
    }

    auto begin() noexcept
    {
        return data.begin();
    }

    auto begin() const noexcept
    {
        return data.begin();
    }

    auto end() noexcept
    {
        return data.end();
    }

    auto end() const noexcept
    {
        return data.end();
    }

    friend bool operator==(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return lhs.data == rhs.data;
    }

    friend bool operator!=(binning_bitvector const & lhs, binning_bitvector const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    auto operator[](size_t const i) noexcept
    {
        assert(i < size());
        return data[i];
    }

    auto operator[](size_t const i) const noexcept
    {
        assert(i < size());
        return data[i];
    }

    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
};

} // namespace raptor::benchmark
