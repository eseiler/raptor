// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <filesystem>

#include <cereal/archives/binary.hpp>

#include <seqan3/contrib/sdsl-lite.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <hibf/config.hpp>
#include <hibf/misc/bit_vector.hpp>

#if 0
#    include <raptor/file_reader.hpp>

// std::cout << "64: " << max_bin_elements(64) << '\n';
// std::cout << "8192: " << max_bin_elements(8192) << '\n';
size_t max_bin_elements(size_t const num_bins)
{
    size_t result{};
    raptor::file_reader<raptor::file_types::sequence> file_reader{seqan3::ungapped{19u}, 23u};

    robin_hood::unordered_flat_set<uint64_t> kmers{};
#    pragma omp parallel for schedule(dynamic) num_threads(32) private(kmers)
    for (size_t i = 0; i < num_bins; ++i)
    {
        kmers.clear();

        std::string path = [&]()
        {
            if (num_bins == 64)
                return std::format("/srv/data/seiler/raptor_2020/data/{}/bins/bin_{:0>2}.fasta", num_bins, i);
            else
                return std::format("/srv/data/seiler/raptor_2020/data/{}/bins/bin_{:0>4}.fasta", num_bins, i);
        }();

        auto insert_it = std::inserter(kmers, kmers.end());
        file_reader.hash_into(path, insert_it);

#    pragma omp critical
        result = std::max(result, kmers.size());
    }
    return result;
}
#endif

static constexpr size_t max_elements_64{21'856'552ULL};
static constexpr size_t max_elements_8192{171'234ULL};

template <typename index_t>
void load_index(index_t & index, std::filesystem::path const & path)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};

    iarchive(index);
}

template <typename index_t>
void store_index(index_t && index, std::filesystem::path const & path)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}

struct old_index
{
    size_t bins{};
    size_t technical_bins{};
    size_t bin_size{};
    size_t hash_shift{};
    size_t bin_words{};
    size_t hash_funs{};
    seqan3::contrib::sdsl::bit_vector data{};

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(bins);
        archive(technical_bins);
        archive(bin_size);
        archive(hash_shift);
        archive(bin_words);
        archive(hash_funs);
        archive(data);
    }

    old_index(std::filesystem::path const & path)
    {
        load_index(*this, path);
    }

    old_index() = default;
    old_index(old_index const &) = default;
    old_index(old_index &&) = default;
    old_index & operator=(old_index const &) = default;
    old_index & operator=(old_index &&) = default;
    ~old_index() = default;
};

struct new_ibf
{
    uint32_t version{1};
    size_t bins{};
    size_t technical_bins{};
    size_t bin_size{};
    size_t hash_shift{};
    size_t bin_words{};
    size_t hash_funs{};
    seqan::hibf::bit_vector data{};
    std::vector<size_t> occupancy{};
    bool track_occupancy{false};

    new_ibf() = default;
    new_ibf(new_ibf const &) = default;
    new_ibf(new_ibf &&) = default;
    new_ibf & operator=(new_ibf const &) = default;
    new_ibf & operator=(new_ibf &&) = default;
    ~new_ibf() = default;

    explicit new_ibf(old_index && old) :
        bins{old.bins},
        technical_bins{old.technical_bins},
        bin_size{old.bin_size},
        hash_shift{old.hash_shift},
        bin_words{old.bin_words},
        hash_funs{old.hash_funs}
    {
        size_t const bit_capacity = old.data.bit_capacity();
        size_t const bit_size = old.data.bit_size();
        size_t const bytes = bit_capacity / 8;

        data.resize(bit_capacity);
        std::memcpy(data.data(), old.data.data(), bytes);
        data.resize(bit_size);

        if (bit_capacity % 64 != 0)
            std::cerr << "Bit capacity is not a multiple of 64.\n";
        if (bit_capacity < bit_size)
            std::cerr << "Bit size is bigger than bit capacity.\n";
        if (bit_size != bins * bin_size)
            std::cerr << "Bit size is not equal to bin size * number of bins.\n";
        if (bit_size != data.size())
            std::cerr << "Bit size is not equal to data size.\n";
        if (bins != technical_bins)
            std::cerr << "Number of bins is not equal to number of technical bins.\n";

        occupancy.resize(bins);
    }

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(version);
        archive(bins);
        archive(technical_bins);
        archive(bin_size);
        archive(hash_shift);
        archive(bin_words);
        archive(hash_funs);
        archive(data);
        archive(occupancy);
        archive(track_occupancy);
    }
};

struct raptor_index
{
    uint32_t version{3u};
    uint64_t window_size{23u};
    seqan3::shape shape{seqan3::ungapped{19u}};
    uint8_t parts{1u};
    std::vector<std::vector<std::string>> bin_path{};
    bool is_hibf{false};
    double fpr{};
    seqan::hibf::config config{};
    new_ibf ibf{};

    template <typename archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(version);
        archive(window_size);
        archive(shape);
        archive(parts);
        archive(bin_path);
        archive(fpr);
        archive(is_hibf);
        archive(config);
        archive(ibf);
    }
};

struct params
{
    std::filesystem::path input_path{};
    std::filesystem::path output_path{};
    size_t bins{};
};

std::vector<std::vector<std::string>> generate_bin_path(size_t const num_bins)
{
    std::vector<std::vector<std::string>> result{};

    for (size_t i = 0; i < num_bins; ++i)
    {
        if (num_bins == 64)
            result.emplace_back(std::vector<std::string>{std::format("bin_{:0>2}.fasta", i)});
        else
            result.emplace_back(std::vector<std::string>{std::format("bin_{:0>4}.fasta", i)});
    }

    return result;
}

double get_fpr(size_t const num_bins, new_ibf const & ibf)
{
    size_t const num_elements = num_bins == 64 ? max_elements_64 : max_elements_8192;
    double const exp_arg = (ibf.hash_funs * num_elements) / static_cast<double>(ibf.bin_size);
    double const log_arg = 1.0 - std::exp(-exp_arg);
    return std::exp(ibf.hash_funs * std::log(log_arg));
}

void run(params params)
{
    new_ibf ibf{old_index{params.input_path}};

    double const fpr = get_fpr(params.bins, ibf);
    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    std::cout << "FPR " << params.bins << ": " << fpr << '\n';
    // b64_4g:   0.006114519521822226 ≈ 0.61%
    // b64_8g:   0.00159148654594952  ≈ 0.16%
    // b8192_4g: 0.006147534720854324 ≈ 0.62%
    // b8192_8g: 0.001600259014483143 ≈ 0.16%

    std::vector<std::vector<std::string>> bin_path = generate_bin_path(params.bins);

    seqan::hibf::config config{.input_fn = [](size_t const, seqan::hibf::insert_iterator) {},
                               .number_of_user_bins = params.bins,
                               .number_of_hash_functions = ibf.hash_funs,
                               .maximum_fpr = fpr};
    config.validate_and_set_defaults();

    if (config.number_of_user_bins != params.bins)
        std::cerr << "number_of_user_bins != params.bins.\n";

    raptor_index index{.bin_path = std::move(bin_path), //
                       .fpr = fpr,
                       .config = config,
                       .ibf = std::move(ibf)};

    store_index(index, params.output_path);
}

int main()
{
    for (auto bins : {64u, 8192u})
    {
        for (auto size : {4u, 8u})
        {
            run({.input_path =
                     std::format("/srv/data/seiler/raptor_2020/indices/ibf_w23_k19_b{}_{}gb.index", bins, size),
                 .output_path =
                     std::format("/srv/data/seiler/raptor_2020/raptor/build/out/ibf_w23_k19_b{}_{}gb.updated.index",
                                 bins,
                                 size),
                 .bins = bins});
        }
    }

    return 0;
}
