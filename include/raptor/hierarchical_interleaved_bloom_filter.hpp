// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------

#pragma once

#include <ranges>
#include <unordered_map>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/unordered_map.hpp>

#ifndef RAPTOR_HIBF_HAS_COUNT
#    define RAPTOR_HIBF_HAS_COUNT 0
#endif

namespace raptor
{

/*!\brief The HIBF binning directory. A data structure that efficiently answers set-membership queries for multiple
 *        bins.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See
 *                           [seqan3::data_layout](https://docs.seqan.de/seqan/3.0.3/group__submodule__dream__index.html#gae9cb143481c46a1774b3cdf5d9fdb518).
 * \see [seqan3::interleaved_bloom_filter][1]
 * \details
 *
 * This class improves the [seqan3::interleaved_bloom_filter][1] by adding additional bookkeeping that allows
 * to establish a hierarchical structure. This structure can then be used to split or merge user bins and distribute
 * them over a variable number of technical bins. In the [seqan3::interleaved_bloom_filter][1], the number of user bins
 * and technical bins is always the same. This causes performance degradation when there are many user bins or the user
 * bins are unevenly distributed.
 *
 * # Terminology
 *
 * ## Technical Bin
 * A Technical Bin represents an actual bin in the binning directory. In the IBF, it stores its kmers in a single Bloom
 * Filter (which is interleaved with all the other BFs).
 *
 * ## User Bin
 * The user may impose a structure on his sequence data in the form of logical groups (e.g. species). When querying the
 * IBF, the user is interested in an answer that differentiates between these groups.
 *
 * # Hierarchical Interleaved Bloom Filter (HIBF)
 *
 * In constrast to the [seqan3::interleaved_bloom_filter][1], the user bins may be split across multiple technical bins
 * , or multiple user bins may be merged into one technical bin. When merging multiple user bins, the HIBF stores
 * another IBF that is built over the user bins constituting the merged bin. This lower-level IBF can then be used
 * to further distinguish between merged bins.
 *
 * In this example, user bin 1 was split into two technical bins. Bins 3, 4, and 5 were merged into a single technical
 * bin, and another IBF was added for the merged bin.
 * \image html hibf.svg width=90%
 *
 * The individual IBFs may have a different number of technical bins and differ in their sizes, allowing an efficient
 * distribution of the user bins.
 *
 * ## Querying
 * To query the Hierarchical Interleaved Bloom Filter for values, call
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent() and use the returned
 * hibf::hierarchical_interleaved_bloom_filter::membership_agent.
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 *
 * To count the occurrences in each user bin of a range of values in the Hierarchical Interleaved Bloom Filter, call
 * hibf::hierarchical_interleaved_bloom_filter::counting_agent() and use
 * the returned hibf::hierarchical_interleaved_bloom_filter::counting_agent_type.
 *
 * ## Thread safety
 *
 * The Interleaved Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * [1]: https://docs.seqan.de/seqan/3.0.3/classseqan3_1_1interleaved__bloom__filter.html
 */
template <seqan3::data_layout data_layout_mode_ = seqan3::data_layout::uncompressed>
class hierarchical_interleaved_bloom_filter
{
public:
    // Forward declaration
    class user_bins;

    // Forward declaration
    class membership_agent;

#if RAPTOR_HIBF_HAS_COUNT
    // Forward declaration
    template <std::integral value_t>
    class counting_agent_type;
#endif

    //!\brief Indicates whether the Interleaved Bloom Filter is compressed.
    static constexpr seqan3::data_layout data_layout_mode = data_layout_mode_;

    //!\brief The type of an individual Bloom filter.
    using ibf_t = seqan3::interleaved_bloom_filter<data_layout_mode_>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hierarchical_interleaved_bloom_filter() = default;                                              //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter const &) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter const &) = default;                        //!< Defaulted.
    hierarchical_interleaved_bloom_filter(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    hierarchical_interleaved_bloom_filter &
    operator=(hierarchical_interleaved_bloom_filter &&) = default; //!< Defaulted.
    ~hierarchical_interleaved_bloom_filter() = default;            //!< Defaulted.

    //!\}

    //!\brief The individual interleaved Bloom filters.
    std::vector<ibf_t> ibf_vector;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the next IBF.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`.
     * If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * If `j != i` is returned, there is a lower level IBF, bin `b` is a merged bin, and `j` is the ID of the lower
     * level IBF in ibf_vector.
     */
    std::vector<std::vector<int64_t>> next_ibf_id;

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the previous IBF.
     * \details
     * Assume we look up an IBF `i`, i.e. `previous_ibf_id[i]`.
     * If `previous_ibf_id[ibf_idx]=ibf_idx` is returned, there is no higher level IBF, and IBF 'i' is thus the root IBF.
     * If `j != number_ibfs` is returned, there is a higher level IBF , with `j` is the ID of the higher
     * level IBF in ibf_vector and bin `b` is the corresponding merged bin in this higher level IBF.
     * \author Myrthe
     */
    std::vector<std::tuple<size_t, size_t>> previous_ibf_id;

    //!\brief The table with FPR or occupancy values, [ibf_idx][bin_idx]
    std::vector<std::vector<double>> fpr_table;

    //!\brief The table with FPR or occupancy values, [ibf_idx][bin_idx]
    std::vector<std::vector<size_t>> occupancy_table; // again it may make more sense to have this as part of ibf_vector, interleaved_bloom_filter does not belong to this project.

    //!\brief Maximum false positive rate for any user bin.
    double fpr_max;

    //!\brief Maximum number of bins of an IBF.
    double t_max;

    //!\brief K-mer size.
    uint8_t k;

    //!\brief minimizer shape
    bool compute_minimiser{false};

    //!\brief The number of hash functions for the IBFs.
    size_t num_hash_functions;

    //!\brief Maximum (or optimal) number of k-mers (bin counts) that IBF could hold at that moment. The vector contains tuples of (ibf_size, ibf_idx), sorted by ibf size.
    std::vector<std::tuple<size_t, size_t>> ibf_sizes;

    //!\brief The underlying user bins.
    user_bins user_bins;

    //!\brief Returns a membership_agent to be used for counting.
    membership_agent membership_agent() const
    {
        return typename hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent{*this};
    }

    /*!\brief Returns the  number of IBFs
     * \author Myrthe Willemsen
     */
    size_t ibf_count(){
        assert(ibf_vector.size() == next_ibf_id.size()); // temporary
        assert(ibf_vector.size() == occupancy_table.size()); // temporary
        assert(ibf_vector.size() == previous_ibf_id.size()); // temporary
        assert(ibf_vector.size() == fpr_table.size()); // temporary
        return ibf_vector.size();
    }

    /*!\brief Returns if a TB is a merged bin.
     * \details If the filename index is smaller than 0, i.e. undefined, it means that the TB is MB.
     * Alternative to current implementation: use next_ibf_id: assume we look up a bin `b` in IBF `i`, i.e. `next_ibf_id[i][b]`. If `i` is returned, there is no lower level IBF, bin `b` is hence not a merged bin.
     * \param[in] ibf_idx
     * \param[in] bin_idx
     * \author Myrthe Willemsen
    */
    bool is_merged_bin(size_t ibf_idx, size_t bin_idx){
        auto const current_filename_index = user_bins.filename_index(ibf_idx, bin_idx);
        assert(ibf_idx < next_ibf_id.size());
        assert(bin_idx < next_ibf_id[ibf_idx].size());
        if (next_ibf_id[ibf_idx][bin_idx] != (int64_t) ibf_idx and next_ibf_id[ibf_idx][bin_idx] != -1){
            assert(current_filename_index < 0);
            return true;
        }else{ return false;}
    }

    //!\brief Updates the tmax based on the number of user bins.
    void update_tmax(){
        auto next_multiple_of_64 = [](size_t value) {return ((value + 63) >> 6) << 6;};   // Takes a number and calculates the next multiple of 64
        size_t user_bin_count = user_bins.size();
        t_max = next_multiple_of_64(std::ceil(std::sqrt(user_bin_count)));
    }


    /*!\brief Returns the number of bins a TB is split into
     * \param[in] ibf_idx
     * \param[in] bin_idx
     * \author Myrthe Willemsen
    */
    size_t is_split_bin(size_t ibf_idx, size_t bin_idx){ // naming is not great
        size_t number_of_bins = 1;
        if (not is_merged_bin(ibf_idx, bin_idx)){ // the function assumes that only leaf bins can be split.
            std::string filename = user_bins[ibf_idx][bin_idx];
            if (filename != ""){number_of_bins = std::get<2>(user_bins.find_filename(filename));}
        }
        return number_of_bins;
    }

    /*!\brief Remove adjacent TBs from the IBF.
     * \param[in] ibf_idx, bin_idx The indices of the TB to be removed.
     * \param[in] number_of_bins. 1 by default, but can be more for split bins.
     * \author Myrthe Willemsen
    */
    void delete_tbs(size_t ibf_idx, size_t bin_idx, size_t number_of_bins=1){
        auto& ibf = ibf_vector[ibf_idx]; // select the IBF
        for (size_t offset=0; offset < number_of_bins; offset++){ // update FPR table and occupancy=#kmer table.
            seqan3::bin_index const bin_index{bin_idx + offset};
            ibf.clear(bin_index);
            fpr_table[ibf_idx][bin_idx+offset] = 0;
            occupancy_table[ibf_idx][bin_idx+offset] = 0;
            next_ibf_id[ibf_idx][bin_idx+offset] = -1; // Try if we can change this to the current ibf rather than -1, or does this itnerfere with the attach subtree function?
        }
    }

    /*!\brief Return the filenames of all user bins contained by a certain IBF, recursively.
     * \param[in] ibf_idx
     * \return a set of strings with filenames
     * \author Myrthe Willemsen
    */
    std::set<std::string> filenames_children(size_t ibf_idx){
        std::set<std::string> filenames;
        for (size_t bin_idx=0; bin_idx < occupancy_table[ibf_idx].size(); ++bin_idx){ // IMPROVEMENT could also be ibf_vector instead of occupancy_table
            if (is_merged_bin(ibf_idx, bin_idx)){ // recursively add filenames for merged bins.
                std::set<std::string> filenames_mb = filenames_children(next_ibf_id[ibf_idx][bin_idx]); // add this set to filenames
                std::merge(filenames.begin(), filenames.end(), filenames_mb.begin(), filenames_mb.end(), std::inserter(filenames, filenames.begin())); // IMPROVEMENT you could use a robin hood set and a lambda function, and make the function recursive.
            }else if(user_bins[ibf_idx][bin_idx] != ""){ // Otherwise an empty string is inserted, which causes downstream problems.
                filenames.insert(user_bins[ibf_idx][bin_idx]);
            }
        }
        return filenames; // set of strings with filenames
    }

    /*!\brief Return the IBF indices of all IBFs in the subtree of a given IBF.
     * \remark the returned child_indices does not include ibf_idx itself.
     * \param[in] ibf_idx the (higher level) IBF of which the lower level IBF indices have to be computed.
     * \param[out] child_indices all (lower level) IBF indices
     * \return child_indices
    * \author Myrthe Willemsen
    */
    void ibf_indices_childeren_(size_t const ibf_idx, std::vector<size_t> & child_indices)
    {
        for (size_t bin_idx=0; bin_idx < next_ibf_id[ibf_idx].size(); ++bin_idx){
            if (is_merged_bin(ibf_idx, bin_idx)){ // recursively for merged bins.
                auto ibf_idx_mb = next_ibf_id[ibf_idx][bin_idx];
                child_indices.push_back(ibf_idx);
                ibf_indices_childeren_(ibf_idx_mb, child_indices);
            }
        }
    }

    std::vector<size_t> ibf_indices_childeren(size_t ibf_idx)
    {
        std::vector<size_t> child_indices{};
        ibf_indices_childeren_(ibf_idx, child_indices);
        return child_indices;
    }


    /*!\brief initialize the previous_ibf_id table
     * \details loop over the (key, item) entries in next_ibf_id, and store the reverse (item, key).
     * \remark If previous_ibf_id[ibf_idx]=ibf_idx, then ibf_idx=root_idx, meaning the IBF at ibf_idx is the root_ibf.
     * \author Myrthe Willemsen
     */
    void initialize_previous_ibf_id(){
        previous_ibf_id.resize(ibf_vector.size());
        std::fill(previous_ibf_id.begin(), previous_ibf_id.end(), std::make_tuple(0,0)); // Initialize with root idx (i.e. 0).
        for (uint64_t ibf_id_high=0; ibf_id_high < next_ibf_id.size(); ++ibf_id_high){
            for (uint64_t bin_idx=0;  bin_idx < next_ibf_id[ibf_id_high].size(); ++bin_idx){
                uint64_t ibf_id_low = next_ibf_id[ibf_id_high][bin_idx]; // -1 in the next_ibf_id will be converted to a very high number, but currently it works out, because ibf_id_low != ibf_id_high still holds.
                if (ibf_id_low != ibf_id_high and ibf_id_low < previous_ibf_id.size()){ // For leaf bins ibf_id_low equals ibf_id_high. For these we should not overwrite previous_ibf_id the entry, because we will lose the information on its parent IBF.
                     previous_ibf_id[ibf_id_low] = std::make_tuple((size_t) ibf_id_high, (size_t) bin_idx);
                }
            }
        }
    }

    /*!\brief Approximates the false positive rate
     * \details Calculates the approximate probability of a false positive hit of a k-mer for a technical bin, based on the bin size, the number of kmers and the number of hash functions.
     * The second function also uses the number of split bins, and returns the FPR for the user bin that was split among those bin. The probability to obtain a false positive answer for the user bin when querying those s techinical bins is 1 - (1 - fpr_tb)^s.
     * The third function also accounts for the .  number_deleted_kmers/(alphabet_size**k) is the probabibility that a too be searched kmer is the same as a deleted kmer (that is still present in the bloom filter)
     * \ref
     * Find the formulas in the corresponding publication
     * \param[in] bin_size size, in number of bits, of the bloom filter.
     * \param[in] kmer_count k-mer count, i.e. the number of kmers stored in the technical bin.
     * \param[in] split split bins: number of bins the among which the n k-mers are divided.
     * \author Myrthe Willemsen
     * \comment It might be computationally faster to work with log2 fpr values, to prevent the power to the h
     */

    double compute_fpr(size_t const bin_size, size_t const kmer_count)
    {
        double const exp_arg = (num_hash_functions * kmer_count) / static_cast<double>(bin_size);
        double const log_arg = 1.0 - std::exp(-exp_arg);
        return std::exp(num_hash_functions * std::log(log_arg));
    }

    double compute_fpr(size_t const bin_size, size_t const kmer_count, size_t const split)
    {
        double const fpr_tb = compute_fpr(bin_size, (kmer_count + split - 1) / split); // fpr for a technical bin. Divide the k-mer count by the number of technical bins
        return 1.0 - std::exp(std::log1p(-fpr_tb) * split);
    }


    /*!\brief Update false positive rate of a TB.
     * \details Calculates the approximate false positive rate of a certain technical bin, given its indices in the HIBF, using the approximate_fpr function and updates .
     * \param[in] ibf_idx index of IBF in HIBF
     * \param[in] bin_idx index of TB in IBF
     * \param[in] number_of_bins index of TB in IBF.
     * \author Myrthe
     */
    double update_fpr(size_t ibf_idx, size_t bin_idx, size_t number_of_bins = 1)
    {   auto kmer_count = occupancy_table[ibf_idx][bin_idx];
        double fpr;
        if (kmer_count==0){
           fpr=0;
        }else{
            auto& ibf = ibf_vector[ibf_idx]; //  select the IBF
            auto bin_size = ibf.bin_size();
            fpr = compute_fpr(bin_size, kmer_count * number_of_bins, // the k-mer count is already over a single bin, so multiply this by the number of bins.
                              number_of_bins); // if a split bin, i think it would be best to store the joint fpr, e.g. of both bins plus multiple testing. We can assume that the occupancy is alomst equal for the tbs among which a UB was split.
        }
        for (size_t offset=0; offset < number_of_bins; ++offset){  //loop over split bins and add multiple testing correction
            fpr_table[ibf_idx][bin_idx+offset] = fpr;
        }
        assert(fpr<1);
        assert(fpr>=0);
        return fpr;
    }

    /*!\brief Obtain the FPR for a technical bin from the FPR table.
     * \author Myrthe Willemsen
    */
    double get_fpr(size_t ibf_idx, size_t const bin_idx){
        return(fpr_table[ibf_idx][bin_idx]);
    }

    /*!\brief Update the occupancy table entry for a specific user bin or adjacent technical bins
     * \param[in] kmer_count the new k-mer count.
     * \param[in] ibf_idx, bin_idx, number_of_bins indices of the technical bin of interest
     * \author Myrthe Willemsen
    */
    void update_occupancy_table(size_t const kmer_count, size_t const ibf_idx, size_t const bin_idx, size_t const number_of_bins = 1){
        if (kmer_count){ // Otherwise ones are placed in the occupancy table
            auto const occupancy_per_bin = std::max<size_t>(1, kmer_count/number_of_bins); // integer devision
            for (size_t offset=0; offset < number_of_bins; offset++){ // update FPR table and occupancy=#kmer table. Perhaps do this before inserting.
                occupancy_table[ibf_idx][bin_idx + offset] += occupancy_per_bin;
            }
        }

    }

    /*!\brief Obtain the number of k-mers that an UB corresponding to the given filename contains.
     * \details The UB might be stored in multiple leaf-bins (LBs). Those occupancies need to be obtained from the occupancy table.
     * \param[in] filename the filename of the UB/sample.
     * \return occupancy_per_file, the number of k-mers that the UB contains.
    * \author Myrthe Willemsen
    */
    size_t get_occupancy_file(std::string const & filename) const
    {
        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = user_bins.find_filename(filename);
        size_t const ibf_idx = std::get<0>(index_triple);
        size_t const bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        size_t occupancy_per_file = 0;
        for (size_t offset=0; offset < number_of_bins; offset++){ // update FPR table and occupancy=#kmer table. Perhaps do this before inserting.
            occupancy_per_file += occupancy_table[ibf_idx][bin_idx + offset];
        }
        return occupancy_per_file;
    }

    /*!\brief Resizes an IBF in the `ibf_vector` and the supporting tables.
     * \details The IBF itself, and the the supporting tables `occupancy_table`, `fpr_table`,
     * `ibf_bin_to_filename_position` are resized according to the new `bin_count`.
     * The latter is a private variable of the user_bins,
     * and is thus resized using a separate function that is part of that class.
     * Note that one does not have to resize the previous_ibf_id.
     * \param[in] ibf_idx index of IBF in HIBF
     * \param[in] bin_count the number of bins to which the IBF has to be resized to.
    * \author Myrthe Willemsen
    */
    void resize_ibf(size_t ibf_idx, size_t bin_count){ //details resize Ibf and datastructures. perhaps move to HIBF.hpp
        ibf_vector[ibf_idx].increase_bin_number_to(seqan3::bin_count{bin_count});
        resize_ibf_occupancy_table(ibf_idx, bin_count);
        resize_ibf_fpr_table(ibf_idx, bin_count);
        user_bins.resize_ibf_filename_positions(ibf_idx, bin_count); // perhaps also update ibf_bin_to_filename_position
        size_t original_size = next_ibf_id[ibf_idx].size();
        next_ibf_id[ibf_idx].resize(bin_count);
        std::fill(next_ibf_id[ibf_idx].begin() + original_size, next_ibf_id[ibf_idx].end(), -1); // initialize the new slots with -1.
    }

    void resize_ibf_occupancy_table(size_t ibf_idx, size_t new_bin_count){
        assert(new_bin_count >= occupancy_table[ibf_idx].size()); // check that new bin count is larger than the current size.
        occupancy_table[ibf_idx].resize(new_bin_count);
    }

    void resize_ibf_fpr_table(size_t ibf_idx, size_t new_bin_count){
        assert(new_bin_count >= fpr_table[ibf_idx].size()); // check that new bin count is larger than the current size.
        fpr_table[ibf_idx].resize(new_bin_count);
    }

    /*!\brief Obtain the number of k-mers that an UB corresponding to the given filename contains.
     * \details The UB might be stored in multiple leaf-bins (LBs). Those occupancies need to be obtained from the occupancy table.
     * \param[in] filename the filename of the UB/sample.
     * \return occupancy_per_file, the number of k-mers that the UB contains.
     * \author Myrthe Willemsen
    */
    double ibf_max_kmers(size_t ibf_idx) const{ // get fpr of a TB in a certain IBF. Makes more sense if this is part of interleaved_bloom_filter.hpp, but that is not part of the project. Myrthe. or bin_index const bin_idx
        auto const& ibf = ibf_vector[ibf_idx]; //  select the IBF       or hibf_ptr->ibf_vector.[ibf_idx] OR  auto& ibf = index.ibf().ibf_vector[ibf_idx]
        size_t const bin_size = ibf.bin_size();
        size_t const hash_funs = ibf.hash_function_count();//arguments.hash
        return approximate_kmer_capacity(bin_size, hash_funs);
    }

    /*!\brief Maximum number of k-mers that fits in a TB
     * \details This can be a value specific to an IBF, which has a given bin size. In addition the bin size in bits can be calculated from:  https://github.com/seqan/raptor/blob/7fe02401bb4f191e2ef4e1454ea9a1c7756816ca/src/build/hibf/bin_size_in_bits.cpp can be used to calculate m from fpr, n and h
     * \param[in] m bin size
     * \param[in] h number of hash functions
     * \return occupancy_per_file, the number of k-mers that the UB contains.
     * \author Myrthe Willemsen
     */
    size_t approximate_kmer_capacity(int m, int h) const
    {
        return -std::floor((m*log(1-pow(fpr_max, (double) 1/h))/h)); // TODO check function .
    }

    /*!\brief Initialize a table with bin sizes of all IBFs.
     * \details The function `ibf_max_kmers` is used to calculate the size of each IBF.
     * which are then stored in a table, such that they do not need to be recalculated.
     * The table of ibf_sizes is used for one of the UB insertion methods of the dynamic HIBF.
     * \author Myrthe Willemsen
     */
    void initialize_ibf_sizes(){
        ibf_sizes.clear();
        for (size_t ibf_idx=0; ibf_idx < ibf_vector.size(); ibf_idx++){
            ibf_sizes.push_back(std::make_tuple(ibf_max_kmers(ibf_idx), ibf_idx));
        }
        std::ranges::sort(ibf_sizes);  // Sort by the first element of tuple
    }

    /*!\brief Calculates the number of TBs needed to store a UB in a specific IBF.
     * \details Given an IBF, with a certain bin size, how many technical bins (TBs) are needed to store a certain
     * number of kmers, considering the multiple testing problem?
     * \param[in] ibf_idx index of IBF in HIBF
     * \param[in] kmer_count the number of k-mers that the UB contains.
     * \return number_of_bins, the number of TBs needed to store the UB.
    * \author Myrthe Willemsen
    */
    size_t number_of_bins(size_t ibf_idx, int kmer_count){
        auto& ibf = ibf_vector[ibf_idx]; // select the IBF
        int bin_size = ibf.bin_size();
        int hash_funs = ibf.hash_function_count();
        // alternatively size_t number_of_bins = (kmer_counts + ibf_max_kmers(ibf_idx) - 1)/ ibf_max_kmers(ibf_idx)
        int number_of_bins = std::ceil(kmer_count / static_cast<double>(ibf_max_kmers(ibf_idx))); // first guess without accounting for multiple testing correction
        while (compute_fpr(bin_size, kmer_count, number_of_bins) > fpr_max){
            number_of_bins++;
        }
        assert(number_of_bins);

        return (size_t) number_of_bins;
    }





#if RAPTOR_HIBF_HAS_COUNT
    /*!\brief Returns a counting_agent_type to be used for counting.
     * \tparam value_t The type to use for the counters; must model std::integral.
     */
    template <std::integral value_t = uint16_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
#endif

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(ibf_vector);
        archive(next_ibf_id);
        archive(user_bins);
        archive(occupancy_table);
        archive(fpr_table);
        archive(previous_ibf_id);
        archive(fpr_max);
        archive(k);
        archive(t_max);
        archive(ibf_sizes);
        archive(num_hash_functions);
    }
    //!\endcond
};

/*!\brief Bookkeeping for user and technical bins.
 */
template <seqan3::data_layout data_layout_mode>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::user_bins
{
public:
    //!\brief Contains filenames of all user bins.
    std::vector<std::string> user_bin_filenames;

    //!\brief Maps filenames to their indices in the list of filenames
    std::unordered_map<std::string, uint64_t> filename_to_idx;

    /*!\brief Stores for each filename ID each of its bins and corresponding IBF in the HIBF
     * \details The `filename_index` is mapped to a triple (ibf_idx, bin_idx, number_of_bins)
     */
    std::vector<std::tuple<uint64_t, uint64_t, uint16_t>> filename_position_to_ibf_bin{};

    /*!\brief Stores for each bin in each IBF of the HIBF the ID of the filename.
     * \details
     * Assume we look up a bin `b` in IBF `i`, i.e. `ibf_bin_to_filename_position[i][b]`.
     * If `-1` is returned, bin `b` is a merged bin, and there is no filename, we need to look into the lower level IBF.
     * Otherwise, the returned value `j` can be used to access the corresponding filename `user_bin_filenames[j]`.
     */
    std::vector<std::vector<int64_t>> ibf_bin_to_filename_position{};

    /*!\brief Initiliazes several of the filename datastructures
     * \details Repopulates the `filename_to_idx` and `filename_position_to_ibf_bin` using the \
     * other datastructures `ibf_bin_to_filename_position` and the `user_bin_filenames` vector
     * \author Myrthe Willemsen
     */
    void initialize_filename_position_to_ibf_bin()
    {
        filename_position_to_ibf_bin.resize(user_bin_filenames.size());
        std::ranges::fill(filename_position_to_ibf_bin, std::make_tuple(0u, 0u, 0u));
        filename_to_idx.clear(); // clear the content of `filename_to_idx`
        for (size_t idx{}; idx < user_bin_filenames.size(); ++idx) // repopulate `filename_to_idx` by traversing over the `user_bin_filenames` vector
        {
            std::string filename = user_bin_filenames[idx];
            if (filename.find(".empty_bin") == std::string::npos) filename_to_idx.emplace(filename, idx); // if it is not an empty bin

        }
        for (size_t ibf_idx{}; ibf_idx < ibf_bin_to_filename_position.size(); ++ibf_idx)
        {
            auto const & ibf_data = ibf_bin_to_filename_position[ibf_idx];

            for (size_t bin_idx{}; bin_idx < ibf_data.size(); ++bin_idx)
            {
                int64_t const filename_position = ibf_data[bin_idx];
                if (filename_position == -1) // merged bin
                    continue;
                assert(filename_position >= 0);
                assert(static_cast<size_t>(filename_position) < filename_position_to_ibf_bin.size());

                auto & ibf_bin = filename_position_to_ibf_bin[filename_position];
                if (std::get<2>(ibf_bin)){ // it is a split bin
                    assert(std::get<0>(ibf_bin) == ibf_idx); // prevent that the same user bin exists in two different places
                    ++std::get<2>(ibf_bin);
                }
                else // as a user bin can take multiple bins, this should consist of multiple bin_idx or multiple tuples
                    ibf_bin = std::make_tuple(ibf_idx, bin_idx, 1u);
            }
        }
    }

     //!\brief resizes an `ibf_filename_positions` vector when resizing an IBF.
    void resize_ibf_filename_positions(size_t ibf_idx, size_t new_bin_count){
        assert(new_bin_count >= ibf_bin_to_filename_position[ibf_idx].size()); // assert that new bin count is larger then the size.
        ibf_bin_to_filename_position[ibf_idx].resize(new_bin_count);
    }

    /*!\brief Returns the bin and IBF indices within the HIBF of a given user bin, specified by filename.
     * \details
     * Should only be used when filename_position_to_ibf_bin has been created, and after checking the filename is present in the filename_to_idx map
     * \author Myrthe Willemsen
     */
    std::tuple<uint64_t, uint64_t, uint16_t> find_filename(std::string const & filename) const
    {
        return filename_position_to_ibf_bin[filename_to_idx.at(filename)];
    }

    /*!\brief Checks if the filename is already present in the HIBF.
     * \author Myrthe Willemsen
     */
    bool exists_filename(const std::string & filename)
    {
    if (filename_to_idx.find(filename) == filename_to_idx.end())
        return false; // filename/user bin is not yet present in HIBF
    else
        return true; // filename/user bin does already exist in HIBF
    }

    /*!\brief update the filename datastructures with a new filename, e.g. when a new user bin has been added
     * \details this function is used for insertions, to update all different user bin datastructure with the new filename.
     * \param [in] filename the filename of the new user bin
     * \param [in] index_triple the ibf_idx, bin_idx and number_of_bins.
     * \author Myrthe Willemsen
     */
    void update_filename_indices(std::string const & filename, std::tuple <uint64_t, uint64_t, uint16_t> & index_triple){
        size_t const ibf_idx = std::get<0>(index_triple);
        size_t const bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        size_t filename_position;
        if (exists_filename(filename)){ // if the user bin already existed, e.g. in the case of sequence insertions, we have to change the existing etnries.
            filename_position = filename_to_idx.at(filename);
            filename_position_to_ibf_bin[filename_position] = index_triple;
            ibf_bin_to_filename_position[ibf_idx][bin_idx] = filename_position ;
        }else{
            filename_position_to_ibf_bin.emplace_back(ibf_idx, bin_idx, number_of_bins);
            filename_position = user_bin_filenames.size();
            filename_to_idx.emplace(filename, filename_position);
            user_bin_filenames.push_back(filename); // make sure this is after the previous methods, such that they refer to the correct idnex.
        }
        for (size_t offset=0; offset < number_of_bins; offset++){
            ibf_bin_to_filename_position[ibf_idx][bin_idx + offset] = filename_position; // do this for all split bins.
        }
    }

    /*!\brief Deletes a filename from the user bin datastructures
     * \details Deletes a filename from `user_bin_filenames`, `filename_to_idx` (in which the filename is set to 'empty_bin',
     * The other two datastructures (`filename_position_to_ibf_bin`, `ibf_bin_to_filename_position`) do not have to be changed now, and will be updated
     * in a future (partial) rebuild operation).
     * \author Myrthe Willemsen
     */
    void delete_filename(const std::string & filename){
        user_bin_filenames[filename_to_idx.at(filename)] ="x.empty_bin"; // One could put the original cardinality or ibf size instead of x, although this does not provide any added value at the moment.
        filename_to_idx.erase(filename);
    }

    //!\brief Returns the number of managed user bins.
    size_t num_user_bins() const noexcept
    {
        return user_bin_filenames.size();
    }

    //!\brief Changes the number of managed IBFs.
    void set_ibf_count(size_t const size)
    {
        ibf_bin_to_filename_position.resize(size);
    }

    //!\brief Changes the number of managed user bins.
    void set_user_bin_count(size_t const size)
    {
        user_bin_filenames.resize(size);
    }

    //!\brief Gets the user bin count.
    size_t size()
    {
        return user_bin_filenames.size();
    }

    //!\brief Returns a vector containing user bin indices for each bin in the `idx`th IBF.
    std::vector<int64_t> & bin_indices_of_ibf(size_t const idx)
    {
        return ibf_bin_to_filename_position[idx];
    }

    //!\brief Returns the filename of the `idx`th user bin.
    std::string & filename_of_user_bin(size_t const idx)
    {
        return user_bin_filenames[idx];
    }

    //!\brief For a pair `(a,b)`, returns a const reference to the filename of the user bin at IBF `a`, bin `b`.
    std::string const & operator[](std::pair<size_t, size_t> const & index_pair) const
    {
        return user_bin_filenames[ibf_bin_to_filename_position[index_pair.first][index_pair.second]];
    }

    /*!\brief Returns a view over the user bin filenames for the `ibf_idx`th IBF.
     *        An empty string is returned for merged bins.
     */
    auto operator[](size_t const ibf_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx]
             | std::views::transform(
                   [this](int64_t i)
                   {
                       if (i == -1)
                           return std::string{};
                       else
                           return user_bin_filenames[i];
                   });
    }

    //!\brief Returns the filename index of the `ibf_idx`th IBF for bin `bin_idx`.
    int64_t filename_index(size_t const ibf_idx, size_t const bin_idx) const
    {
        return ibf_bin_to_filename_position[ibf_idx][bin_idx];
    }

    /*!\brief Writes all filenames to a stream. Index and filename are tab-separated.
     * \details
     * 0	\<path_to_user_bin_0\>
     * 1	\<path_to_user_bin_1\>
     */
    template <typename stream_t>
    void write_filenames(stream_t & out_stream) const
    {   size_t position{};
        std::string line{};
        for (auto const & filename : user_bin_filenames)
        {
            line.clear();
            line = '#';
            line += std::to_string(position);
            line += '\t';
            line += filename;
            line += '\n';
            out_stream << line;
            ++position;
        }
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <typename archive_t>
    void serialize(archive_t & archive)
    {
        archive(user_bin_filenames);
        archive(ibf_bin_to_filename_position);
        archive(filename_position_to_ibf_bin);
        archive(filename_to_idx);
    }
    //!\endcond
    };

/*!\brief Manages membership queries for the hibf::hierarchical_interleaved_bloom_filter.
 * \see hibf::hierarchical_interleaved_bloom_filter::user_bins::filename_of_user_bin
 * \details
 * In contrast to the [seqan3::interleaved_bloom_filter][1], the result will consist of indices of user bins.
 */
template <seqan3::data_layout data_layout_mode>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::membership_agent
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive membership querying.
    template <std::ranges::forward_range value_range_t>
    void bulk_contains_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<uint16_t>();
        auto & result = agent.bulk_count(values);

        uint16_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];

            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold) // if ibf contains the query.
                    bulk_contains_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold); //recursive
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer.emplace_back(current_filename_index);
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    membership_agent() = default;                                     //!< Defaulted.
    membership_agent(membership_agent const &) = default;             //!< Defaulted.
    membership_agent & operator=(membership_agent const &) = default; //!< Defaulted.
    membership_agent(membership_agent &&) = default;                  //!< Defaulted.
    membership_agent & operator=(membership_agent &&) = default;      //!< Defaulted.
    ~membership_agent() = default;                                    //!< Defaulted.

    /*!\brief Construct a membership_agent for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit membership_agent(hibf_t const & hibf) : hibf_ptr(std::addressof(hibf))
    {}
    //!\}

    //!\brief Stores the result of bulk_contains().
    std::vector<int64_t> result_buffer;

    /*!\name Lookup
     * \{
     */
    /*!\brief Determines set membership of given values, and returns the user bin indices of occurrences.
     * \param[in] values The values to process; must model std::ranges::forward_range.
     * \param[in] threshold Report a user bin if there are at least this many hits.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::membership_agent for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values, size_t const threshold) & noexcept
    {
        assert(hibf_ptr != nullptr);

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        result_buffer.clear();

        bulk_contains_impl(values, 0, threshold);

        std::ranges::sort(result_buffer);

        return result_buffer;
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] std::vector<int64_t> const & bulk_contains(value_range_t && values,
                                                             size_t const threshold) && noexcept = delete;
    //!\}
};


#if RAPTOR_HIBF_HAS_COUNT
/*!\brief Manages counting ranges of values for the hibf::hierarchical_interleaved_bloom_filter.
 */
template <seqan3::data_layout data_layout_mode>
template <std::integral value_t>
class hierarchical_interleaved_bloom_filter<data_layout_mode>::counting_agent_type
{
private:
    //!\brief The type of the augmented hierarchical_interleaved_bloom_filter.
    using hibf_t = hierarchical_interleaved_bloom_filter<data_layout_mode>;

    //!\brief A pointer to the augmented hierarchical_interleaved_bloom_filter.
    hibf_t const * const hibf_ptr{nullptr};

    //!\brief Helper for recursive bulk counting.
    template <std::ranges::forward_range value_range_t>
    void bulk_count_impl(value_range_t && values, int64_t const ibf_idx, size_t const threshold)
    {
        auto agent = hibf_ptr->ibf_vector[ibf_idx].template counting_agent<value_t>();
        auto & result = agent.bulk_count(values);

        value_t sum{};

        for (size_t bin{}; bin < result.size(); ++bin)
        {
            sum += result[bin];
            auto const current_filename_index = hibf_ptr->user_bins.filename_index(ibf_idx, bin);

            if (current_filename_index < 0) // merged bin
            {
                if (sum >= threshold)
                    bulk_count_impl(values, hibf_ptr->next_ibf_id[ibf_idx][bin], threshold);
                sum = 0u;
            }
            else if (bin + 1u == result.size() ||                                                    // last bin
                     current_filename_index != hibf_ptr->user_bins.filename_index(ibf_idx, bin + 1)) // end of split bin
            {
                if (sum >= threshold)
                    result_buffer[current_filename_index] = sum;
                sum = 0u;
            }
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default;                                        //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default;             //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default;                  //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default;      //!< Defaulted.
    ~counting_agent_type() = default;                                       //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing hierarchical_interleaved_bloom_filter.
     * \private
     * \param hibf The hierarchical_interleaved_bloom_filter.
     */
    explicit counting_agent_type(hibf_t const & hibf) :
        hibf_ptr(std::addressof(hibf)),
        result_buffer(hibf_ptr->user_bins.num_user_bins())
    {}
    //!\}

    //!\brief Stores the result of bulk_count().
    seqan3::counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences in each bin for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::forward_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     * \param[in] threshold Do not recurse into merged bins with less than this many hits. Default: 1.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &`, to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a
     * hibf::hierarchical_interleaved_bloom_filter::counting_agent_type for each thread.
     */
    template <std::ranges::forward_range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) & noexcept
    {
        assert(hibf_ptr != nullptr);
        assert(threshold > 0u);
        assert(result_buffer.size() == hibf_ptr->user_bins.num_user_bins());

        static_assert(std::ranges::forward_range<value_range_t>, "The values must model forward_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        std::ranges::fill(result_buffer, static_cast<value_t>(0u));

        bulk_count_impl(values, 0, threshold);

        return result_buffer;
    }

    // `bulk_count` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    template <std::ranges::range value_range_t>
    [[nodiscard]] seqan3::counting_vector<value_t> const & bulk_count(value_range_t && values,
                                                                      size_t const threshold = 1u) && noexcept = delete;
    //!\}
};
#endif // RAPTOR_HIBF_HAS_COUNT

} // namespace raptor
