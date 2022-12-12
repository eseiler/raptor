// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#include <raptor/test/benchmark/binning_bitvector.hpp>

namespace raptor::benchmark
{

template <std::integral value_t>
class counting_vector : public std::vector<value_t>
{
    using base_t = std::vector<value_t>;

public:
    counting_vector() = default;
    counting_vector(counting_vector const &) = default;
    counting_vector & operator=(counting_vector const &) = default;
    counting_vector(counting_vector &&) = default;
    counting_vector & operator=(counting_vector &&) = default;
    ~counting_vector() = default;

    using base_t::base_t;

    template <typename binning_bitvector_t>
    counting_vector & operator+=(binning_bitvector_t const & binning_bitvector)
    {
        for_each_set_bin(binning_bitvector,
                         [this](size_t const bin)
                         {
                             ++(*this)[bin];
                         });
        return *this;
    }

private:
    template <typename binning_bitvector_t, typename on_bin_fn_t>
    void for_each_set_bin(binning_bitvector_t && binning_bitvector, on_bin_fn_t && on_bin_fn)
    {
        assert(this->size() >= binning_bitvector.size());

        auto jump_to_next_1bit = [](size_t & x)
        {
            auto const zeros = std::countr_zero(x);
            x >>= zeros; // skip number of zeros
            return zeros;
        };

        for (size_t bit_pos = 0; bit_pos < binning_bitvector.size(); bit_pos += 64)
        {
            size_t bit_sequence = binning_bitvector.raw_data().get_int(bit_pos);

            for (size_t bin = bit_pos; bit_sequence != 0u; ++bin, bit_sequence >>= 1)
            {
                bin += jump_to_next_1bit(bit_sequence);

                on_bin_fn(bin);
            }
        }
    }
};

template <std::integral value_t>
class counting_vector_invert : public std::vector<value_t>
{
    using base_t = std::vector<value_t>;

public:
    counting_vector_invert() = default;
    counting_vector_invert(counting_vector_invert const &) = default;
    counting_vector_invert & operator=(counting_vector_invert const &) = default;
    counting_vector_invert(counting_vector_invert &&) = default;
    counting_vector_invert & operator=(counting_vector_invert &&) = default;
    ~counting_vector_invert() = default;

    using base_t::base_t;

    template <typename binning_bitvector_t>
    counting_vector_invert & operator+=(binning_bitvector_t const & binning_bitvector)
    {
        for_each_set_bin(binning_bitvector,
                         [this](size_t const bin)
                         {
                             ++(*this)[bin];
                         });
        return *this;
    }

private:
    template <typename binning_bitvector_t, typename on_bin_fn_t>
    void for_each_set_bin(binning_bitvector_t && binning_bitvector, on_bin_fn_t && on_bin_fn)
    {
        assert(this->size() >= binning_bitvector.size());

        auto jump_to_next_1bit = [](size_t & x)
        {
            auto const zeros = std::countr_zero(x);
            x >>= zeros; // skip number of zeros
            return zeros;
        };

        for (size_t bit_pos = 0; bit_pos < binning_bitvector.size(); bit_pos += 64)
        {
            size_t bit_sequence = ~(binning_bitvector.raw_data().get_int(bit_pos));

            if (bit_sequence == 0u)
            {
                for (size_t i = 0; i < 64; ++i)
                    on_bin_fn(i);
            }
            else
            {
                for (size_t bin = bit_pos; bin < bit_pos + 64; ++bin, bit_sequence >>= 1)
                {
                    size_t const next_unset = jump_to_next_1bit(bit_sequence);

                    for (size_t i = 0; i < next_unset; ++i)
                        on_bin_fn(bin++);
                }
            }
        }
    }
};

template <std::integral value_t>
class counting_vector_invert_all : public std::vector<value_t>
{
    using base_t = std::vector<value_t>;

public:
    counting_vector_invert_all() = default;
    counting_vector_invert_all(counting_vector_invert_all const &) = default;
    counting_vector_invert_all & operator=(counting_vector_invert_all const &) = default;
    counting_vector_invert_all(counting_vector_invert_all &&) = default;
    counting_vector_invert_all & operator=(counting_vector_invert_all &&) = default;
    ~counting_vector_invert_all() = default;

    using base_t::base_t;

    template <typename binning_bitvector_t>
    counting_vector_invert_all & operator+=(binning_bitvector_t const & binning_bitvector)
    {
        for_each_set_bin(binning_bitvector,
                         [this](size_t const bin)
                         {
                             --(*this)[bin];
                         });
        return *this;
    }

private:
    template <typename binning_bitvector_t, typename on_bin_fn_t>
    void for_each_set_bin(binning_bitvector_t && binning_bitvector, on_bin_fn_t && on_bin_fn)
    {
        assert(this->size() >= binning_bitvector.size());

        std::ranges::for_each(*this,
                              [](value_t & elem)
                              {
                                  ++elem;
                              });

        auto jump_to_next_1bit = [](size_t & x)
        {
            auto const zeros = std::countr_zero(x);
            x >>= zeros; // skip number of zeros
            return zeros;
        };

        for (size_t bit_pos = 0; bit_pos < binning_bitvector.size(); bit_pos += 64)
        {
            size_t bit_sequence = ~(binning_bitvector.raw_data().get_int(bit_pos));

            for (size_t bin = bit_pos; bit_sequence != 0u; ++bin, bit_sequence >>= 1)
            {
                bin += jump_to_next_1bit(bit_sequence);

                on_bin_fn(bin);
            }
        }
    }
};

} // namespace raptor::benchmark
