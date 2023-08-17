//#include <raptor/search/search_single.hpp> // to make sure index_structure::hibf_compresse exists here.
#include <lemon/list_graph.h> /// Must be first include.

#include <chopper/layout/execute.hpp>
#include <chopper/set_up_parser.hpp>
#include <raptor/update/load_hibf.hpp>
#include <raptor/update/insertions.hpp>
#include <raptor/build/store_index.hpp>
#include <raptor/update/rebuild.hpp>
#include <raptor/build/hibf/chopper_build.hpp>
#include <raptor/build/hibf/create_ibfs_from_chopper_pack.hpp> // when adding this, it recognizes std as raptor::std
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <chopper/data_store.hpp>
#include <chopper/layout/insert_empty_bins.hpp>
#include <chopper/next_multiple_of_64.hpp>
#include <random>
#include "chopper/sketch/execute.hpp"


namespace raptor
{

//!\brief helper function to create random filenames for temporary folders to save the layout in. This is to prevent multiple runs writing to the same file.
std::string generate_random_string(int length) {
    std::string characters = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, characters.length() - 1);

    std::string result;
    for (int i = 0; i < length; ++i) {
        result += characters[dis(gen)];
    }

    return result;
}
const extern std::string tmp_folder = "tmp" + generate_random_string(10); // create a random tmp folder, s.t. multiple runs will not be in each others way.

//!\brief helper function to convert vector of strings to set
std::set<std::string> convert_to_set(std::vector<std::string> vector)
{
    std::set<std::string> set;
    for (std::string string : vector)  if (string != " " and string != "") set.insert(string);
    return set;
}

//!\brief Rebuilds the complete index. Otherwise works similar to partial rebuild.
void full_rebuild(raptor_index<index_structure::hibf> & index,
                  update_arguments const & update_arguments) {
    std::cout << "Full rebuild" << std::flush;
    index.ibf().initialize_ibf_sizes(); // I think this is needed here to prevent assertion statements, but I am not sure.
    //0) Create layout arguments
    chopper::configuration layout_arguments = layout_config(index, update_arguments); // create the arguments to run the layout algorithm with.
    //1) Obtain kmer counts together with the filenames
    index.ibf().ibf_vector = {}; // index.ibf().ibf_vector = std::vector<ibf_type>{}; //remove the original IBF to reduce RAM usage.
    std::vector<std::string> filenames = index.ibf().user_bins.user_bin_filenames; // all filenames
    auto kmer_counts_filenames = get_kmer_counts(index, convert_to_set(filenames));
    //2) call chopper layout on the stored filenames.
    call_layout(kmer_counts_filenames, layout_arguments);
    //3) call hierarchical build.
    raptor_index<index_structure::hibf> new_index{}; //create the empty HIBF of the subtree.
    build_arguments build_arguments = build_config(index, layout_arguments); // create the arguments to run the build algorithm with.
    call_build(build_arguments, new_index, true); // The additional datastructures are therein also created
    index.ibf() = new_index.ibf();
    std::filesystem::remove_all(tmp_folder);
}


/*!\brief Checks if the bin count of an IBF exceeds t_max, and requires a rebuild.
 * \details Triggers a rebuild directly from the function.
 * \param[in] ibf_idx the location of the subtree in the HIBF
 * \param[in] index the original HIBF
 * \param[in] update_arguments the arguments that were passed with the update that was to be done on the HIBF.
 * \return if a rebuild was performed.
 * \author Myrthe
 */
bool check_tmax_rebuild(raptor_index<index_structure::hibf> & index, size_t ibf_idx,
                        update_arguments const & update_arguments){
    if (index.ibf().ibf_vector[ibf_idx].bin_count() > index.ibf().t_max){ // If an ibf grows out of the tolerated t_max, a full rebuild is triggered
        index.ibf().update_tmax(); // first update the t_max
        if (index.ibf().ibf_vector[ibf_idx].bin_count() > index.ibf().t_max){
            if (ibf_idx > 0 and update_arguments.tmax_condition){
                partial_rebuild(index.ibf().previous_ibf_id[ibf_idx], index, update_arguments); // if update_arguments.tmax_condition is false, then the t_max only holds for the top level.
                return true;
            }
            else {
                full_rebuild(index, update_arguments);
                return true;
            }
        }
    }
    return false;
}

/*!\brief Does a partial rebuild for the subtree at index ibf_idx
 * \details The algorithm rebuilds a subtree of the hierarchical interleaved bloom filter completely.
 * (1) It computes the set of filenames belonging to all user bins in the subtree, and stores those in a text file
 * (2) It computes a good layout by calling Choppers's layout algorithm.
 * (3) It rebuilds the subtree by calling hierachical build function.
 * (4) It merges the index of the HIBF of the newly obtained subtree with the original index.
 * \param[in] ibf_idx the location of the subtree in the HIBF
 * \param[in] index the original HIBF
 * \param[in] update_arguments the arguments that were passed with the update that was to be done on the HIBF.
 * \param[in] number_of_splits In how many seperate IBFs must the IBF be split, 2 by default.
 * \author Myrthe
 */
void partial_rebuild(std::tuple<size_t,size_t> index_tuple,
                  raptor_index<index_structure::hibf> & index,
                  update_arguments const & update_arguments, //TODO later add argument with filename + kmer count to be added.
                  int number_of_splits) // instead of setting a default parameter here, add this to the update arguments
{
    std::cout << "Partial rebuild" << std::flush;
    size_t ibf_idx = std::get<0>(index_tuple);
    size_t bin_idx = std::get<1>(index_tuple);
    size_t next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][bin_idx];
    //1.1) Obtain filenames from all lower bins and kmer counts. Perhaps using occupancy table.
    auto filenames_subtree = index.ibf().filenames_children(next_ibf_idx);
    //1.2) Obtain the kmer counts
    auto kmer_counts_filenames = get_kmer_counts(index, filenames_subtree);
    //1.3) Define how to split te IBF, in terms of merged bin indexes and user bin filenames
    auto split_files = split_filenames(kmer_counts_filenames, number_of_splits);
    auto split_idxs = split_mb(index_tuple, index, update_arguments, number_of_splits); // returns the new indices to attach the subtree

    //1.5) Check the IBF bin count exceeds tmax, which might be possible if it has been resized. If tmax is exceeded, another build will be triggered inside the 'check tmax' function
    if (check_tmax_rebuild(index, ibf_idx, update_arguments)) return;

    //1.4) remove IBFs of the to-be-rebuild subtree from the  index.
    remove_ibfs(index, next_ibf_idx);

    for (int split = 0; split < number_of_splits; split++){ // for each subindex that is rebuild.
        if (split_files[split].size()) {
            if (split_files[split].size() == 1){ // if the subtree consists of a single user bin, than the user bin will be directly placed on the higher-level IBF.
                robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
                std::vector<std::string> filename = {std::get<1>(split_files[split][0])};
                raptor::hibf::compute_kmers(kmers, update_arguments, filename);
                insert_into_ibf(kmers, std::make_tuple(ibf_idx, split_idxs[split], 1), index);       // fills the empty bin. no need to insert into parents.
            }else{
                //0) Create layout arguments
                chopper::configuration layout_arguments = layout_config(index, update_arguments, std::to_string(split)); // create the arguments to run the layout algorithm with.
                //2) call chopper layout on the stored filenames.
                call_layout(split_files[split], layout_arguments);
                //3) call hierarchical build.
                raptor_index<index_structure::hibf> subindex{}; //create the empty HIBF of the subtree.
                build_arguments build_arguments = build_config(index, layout_arguments); // create the arguments to run the build algorithm with.
                robin_hood::unordered_flat_set<size_t> root_kmers = call_build(build_arguments, subindex, false); // the last argument sets the 'is_root' parameter to false such that kmers are added to `parent_kmers`. However, this causes that no shuffling takes place reducing the efficiency a bit. This could be further optimized perhaps by adding an extra parameter.
                insert_into_ibf(root_kmers, std::make_tuple(ibf_idx, split_idxs[split], 1), index);       // also fills the (new) MB.
                //5) merge the index of the HIBF of the newly obtained subtree with the original index.
                attach_subindex(index, subindex, std::make_tuple(ibf_idx, split_idxs[split]));
            }
        }
    }
    // update datastructures
    index.ibf().initialize_ibf_sizes();
    index.ibf().user_bins.initialize_filename_position_to_ibf_bin(); // this also updates the filename_to_idx datastructure
    // remove temporary files
    std::filesystem::remove_all(tmp_folder);
}


/*!\brief splits a merged bin
 * \details splits a merged bin into 'number_of_splits', 2 by default, by finding empty bins on the IBF.
 * \param[in] index_tuple the ibf_idx and bin_idx of the merged bin that has reached the FPR limit.
 * \param[in] index the original HIBF
 * \param[in] number_of_splits In how many seperate IBFs must the IBF be split, 2 by default. This may be 1 if it is not split.
 * \author Myrthe Willemsen
 */
std::vector<uint64_t> split_mb(std::tuple<size_t,size_t> index_tuple,
               raptor_index<index_structure::hibf> & index,
               update_arguments const & update_arguments,
               int number_of_splits) // should be size_t
{
    std::vector<uint64_t> tb_idxs(number_of_splits);
    index.ibf().delete_tbs(std::get<0>(index_tuple), std::get<1>(index_tuple));    // Empty the merged bin.
    auto & ibf_idx = std::get<0>(index_tuple);
    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level IBF to serve as new merged bins.
        tb_idxs[split] = find_empty_bin_idx(index, ibf_idx, update_arguments);       // find an empty bin for a new MB on this IBF or resize. Import from insertions.
        size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
        if (tb_idxs[split] == ibf_bin_count){ // current solution when tmax is reached, because there are not sufficient empty bins
                size_t new_ibf_bin_count = new_bin_count(1u, update_arguments.empty_bin_percentage, ibf_bin_count);
                std::cout << "Resize the IBF at index: " << ibf_idx << "\n" << std::flush;
                index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
            }
        else{
            assert(std::get<0>(index_tuple) < index.ibf().occupancy_table.size());
            assert(tb_idxs[split] < index.ibf().occupancy_table[ibf_idx].size());
            index.ibf().occupancy_table[ibf_idx][tb_idxs[split]] = 1000; // set some value to the occupancy table, such that not the same value is found twice
        }
    }
    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level IBF to serve as new merged bins.
        index.ibf().occupancy_table[ibf_idx][tb_idxs[split]] = 0; // set back to 0
    }
    return tb_idxs;
}


/*!\brief Splits the set of user bins in the subtree in `n` similarly sized subsets.
 * \details User bins are devided among `n` subindexes, `n` being th number of splits.
 * The IBF is split on a user bin level, to make sure that no user bin enters both IBFs, which is otherwise difficult because of split bins.
 * Alternatively, one can use the counts for these filenames by calling Chopper count.
 * To take into account sequence similarity, the splitting should be done after rearranging user bins as part of Chopper layout.
 * \param[in] kmer_counts_filenames a vector of tuples with k-mer counts (1) for each filename (2) present in the subtree of the IBF.
 * \param[in] number_of_splits In how many seperate IBFs must the IBF be split, 2 by default.
 * \author Myrthe Willemsen
 */
std::vector<std::vector<std::tuple<size_t, std::string>>> split_filenames(
        std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
        int number_of_splits){
    assert(number_of_splits);
    size_t sum_kmer_count = 0;
    for (auto& tuple : kmer_counts_filenames) sum_kmer_count += std::get<0>(tuple);
    size_t percentile = sum_kmer_count/number_of_splits; // sum over k-mer counts and find when the threshold of the sum is first exceeded.
    size_t cumulative_sum = 0; size_t filename_idx = 0; size_t split_idx =0;
    std::vector<std::vector<std::tuple<size_t, std::string>>> split_filenames;
    for (int split = 0; split < number_of_splits; split++){             // get indices of the empty bins on the higher level to serve as new merged bins.
        while (cumulative_sum < percentile * (split + 1) // while the threshold for the next split not yet exceeded
        && filename_idx < kmer_counts_filenames.size()){ // and while we did not reach the last file
            cumulative_sum += std::get<0>(kmer_counts_filenames[filename_idx]);
            filename_idx += 1;
        }
        if (split==number_of_splits-1) filename_idx = kmer_counts_filenames.size();
        split_filenames.push_back(std::vector(        // create a new vector with filename indices.
            std::ranges::next(kmer_counts_filenames.begin(), split_idx, kmer_counts_filenames.end()),  // std::ranges::next(iterator, number, bound) is the same as iterator + number, but bound: it cannot go out of range.
            std::ranges::next(kmer_counts_filenames.begin(), filename_idx , kmer_counts_filenames.end())));
        split_idx = filename_idx ;
    }
    return split_filenames;
}


/*!\brief Store kmer or minimizer counts for each specified file stored within the HIBF.
 * \details The algorithm obtains the total kmer count for each file and stores those counts together with the
 * filenames to a text file. Output is sorted by count. Alternatively, one can create/extract counts for these filenames by calling Chopper count
 * \param[in] filenames a set of filenames of which the counts should be obtained
 * \param[in] index the original HIBF
 * \return[out] kmer_counts a vector of tuples with kmer counts and filenames.
 * \author Myrthe Willemsen
 */
std::vector<std::tuple<size_t, std::string>> get_kmer_counts(raptor_index<index_structure::hibf> const & index,
                                                             std::set<std::string> const & filenames)
{
    std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames{};
    for (std::string const & filename : filenames)
    {
        if (std::filesystem::path filename_as_path{filename};
        filename_as_path.extension() != ".empty_bin" and filename != "") // construct filename_as_path such that it can be used for comparison without going out of scope
        {
            int const kmer_count = index.ibf().get_occupancy_file(filename);
            kmer_counts_filenames.push_back(std::make_tuple(kmer_count, filename));
    	    std::cout << "filename, kmer_count: " << filename << ", " << kmer_count <<std::endl;
    	}
    }
    std::ranges::sort(kmer_counts_filenames); //  std::ranges::sort(kmer_counts_filenames, std::ranges::greater);
    std::reverse(kmer_counts_filenames.begin(), kmer_counts_filenames.end());
    return kmer_counts_filenames; // array of tuples with filename and k-mer count, sorted descending by kmer count.
}



/*!\brief Creates a configuration object which is passed to chopper's layout algorithm.
 * \param[in] update_arguments the file containing all paths to the user bins for which a layout should be computed
 * \param[in] index the original HIBF
 * \param[in] file_indicator an identifier to distinguish files.
 * \author Myrthe Willemsen
 */
chopper::configuration layout_config(raptor_index<index_structure::hibf> & index,
                                     update_arguments const & update_arguments,
                                     std::string file_indicator){
    std::filesystem::remove_all(tmp_folder);
    std::filesystem::create_directory(tmp_folder);
    chopper::configuration config{};
    //config.data_file = "tmp/temporary_layout" + file_indicator;
    config.output_filename = tmp_folder +"/temporary_layout" + file_indicator + ".txt"; // seqan tmp folder can be used later.
    config.data_file = tmp_folder + "/subtree_bin_paths.txt"; // input of bins paths
    config.sketch_directory = update_arguments.sketch_directory;
    config.disable_rearrangement = not update_arguments.similarity; // indicates whether updates should account for user bin's similarities. This also determines "estimate union"
    config.disable_estimate_union = not update_arguments.similarity;
    config.update_ubs = update_arguments.empty_bin_percentage; // percentage of empty bins drawn from distributions //makes sure to use empty bins.

    index.ibf().update_tmax();
    config.tmax = index.ibf().t_max;
    config.determine_best_tmax = false;

    config.false_positive_rate = index.ibf().fpr_max;
    if (update_arguments.insert_sequences)
        config.false_positive_rate = index.ibf().fpr_max * (1 - update_arguments.empty_bin_percentage); // create an FPR buffer when working with sequence insertions. We cant do this based on 10% extra sequence insertions, you cannot solve that analytically. Here i supppose empty bin perscentage <1

    config.k = index.ibf().k;
    config.num_hash_functions = index.ibf().num_hash_functions;

    config.sketch_bits = update_arguments.sketch_bits;
    config.threads = update_arguments.threads;

    return config;
}


/*!\brief Calls the layout algorithm from the chopper library
* \param[in] layout_arguments configuration object with parameters required for calling the layout algorithm
* \warning an extra enter/return at the start of the bin file will cause segmentation faults in chopper.
* \author Myrthe Willemsen
*/
void call_layout(std::vector<std::tuple<size_t, std::string>> kmer_counts_filenames,
                 chopper::configuration & config){ // layout_arguments
    int exit_code{};

    chopper::layout::layout hibf_layout{};
    std::vector<std::string> filenames;
    std::vector<size_t> kmer_counts;
    for (const auto &entry: kmer_counts_filenames) { // kmer_counts_filenames is sorted decending.
        kmer_counts.push_back(std::get<0>(entry));
        filenames.push_back(std::get<1>(entry));
    }
    std::vector<bool> empty_bins; // A bitvector indicating whether a bin is empty (1) or not (0).
    std::vector<size_t> empty_bin_cum_sizes; // The cumulative k-mer count for the first empty bin until empty bin i
    std::vector<chopper::sketch::hyperloglog> sketches{};

    try {
        try{
            chopper::sketch::toolbox::read_hll_files_into(config.sketch_directory, filenames, sketches); // TODO add 'config' as last argument to  make sure that, if for any sketch, it is not found, that it is created.
        } catch (...) { //temporary alternative solution until chopper code is re-ran.
            std::vector<chopper::sketch::hyperloglog> sketches{};
            chopper::sketch::execute(config, filenames, sketches);
        }

        chopper::layout::insert_empty_bins(empty_bins, empty_bin_cum_sizes,
                          kmer_counts, sketches, filenames, config);

        chopper::data_store store{.false_positive_rate = config.false_positive_rate,
                .hibf_layout = &hibf_layout,
                .kmer_counts = kmer_counts,
                .sketches = sketches,
                .empty_bins = empty_bins,
                .empty_bin_cum_sizes = empty_bin_cum_sizes};

        exit_code |= chopper::layout::execute(config, filenames, store);
    }
    catch (sharg::parser_error const &ext) {    // GCOVR_EXCL_START
        std::cerr << "[CHOPPER ERROR] " << ext.what() << '\n';
    }

}

/*!\brief Creates a configuration object which is passed to hierarchical build function.
* \author Myrthe Willemsen
*/
build_arguments build_config(raptor_index<index_structure::hibf> & index,
                             chopper::configuration layout_arguments){
    build_arguments build_arguments{};
    build_arguments.shape = index.shape(); // shape is more important then k-mer size.
    build_arguments.kmer_size = index.ibf().k;
    build_arguments.window_size =  index.window_size();
    build_arguments.fpr = layout_arguments.false_positive_rate;
    build_arguments.is_hibf = true;
    build_arguments.bin_file = layout_arguments.output_filename;
    build_arguments.threads = layout_arguments.threads;
    return build_arguments;
}

/*!\brief Calls the hierarchical build algorithm
* \param[in] build_arguments configuration object with parameters required for calling the building algorithm
* \param[out] index the newly created index
* \return root_kmers the total set of k-mers stored at the root, and that can be used to store in a merged bin if this regards a subindex.
* \author Myrthe Willemsen
*/
template <seqan3::data_layout data_layout_mode>
robin_hood::unordered_flat_set<size_t> call_build(build_arguments & arguments,
                raptor_index<hierarchical_interleaved_bloom_filter<data_layout_mode>> & index,
                bool is_root){
    hibf::build_data<data_layout_mode> data{};
    robin_hood::unordered_flat_set<size_t> root_kmers = raptor::hibf::create_ibfs_from_chopper_pack(data, arguments, is_root);
    std::vector<std::vector<std::string>> bin_path{};
    for (size_t i{0}; i < data.hibf.user_bins.num_user_bins(); ++i)
        bin_path.push_back(std::vector<std::string>{data.hibf.user_bins.filename_of_user_bin(i)});
    index.ibf() = std::move(data.hibf); //instead of creating the index object here.
    return root_kmers;
}
/*!\brief Helper function that removes the given indices from a vector.
* \author Myrthe Willemsen
*/
template <typename T> void remove_indices(std::vector<size_t> indices_to_remove, std::vector<T> & vector) {
    if (indices_to_remove.size() > 1) std::sort(indices_to_remove.rbegin(), indices_to_remove.rend());   // sort descending to prevent wrong indices being removed on the way. Note that this has resulted in an error in the sdsl library. If it still does, use normal sort and a reverse loop
        for (int i : indices_to_remove) {
            vector.erase(vector.begin() + i);
        }
    }

/*!\brief Prunes subtree from the original HIBF
 * \details One should remove the IBFs in the original index which were part of the subtree that had to be rebuild.
 * `indices_map` is used to map IBF indices of the original HIBF to those in the new HIBF. If the IBF at an old index is not in the new HIBF, it will return -1.
 * If using some sort of splitting, then removing only needs to happen once since both new subindexes share the same original ibfs.
 * \param[in|out] index the original HIBF
 * \param[in] ibf_idx the index of the IBF where the subtree needs to be removed, including the ibf_idx itself.
 * \author Myrthe Willemsen
 */

void remove_ibfs(raptor_index<index_structure::hibf> & index, size_t ibf_idx){
    // Store which original indices in the IBF were the subindex that had to be rebuild?
    std::vector<size_t> indices_to_remove = index.ibf().ibf_indices_childeren(ibf_idx);  // Create a map that maps remaining IBF indices of the original HIBF to their new indices
    indices_to_remove.push_back(ibf_idx);
    std::vector<int> indices_map(index.ibf().ibf_vector.size(), -1);  // The map has the size in number of IBF's of the original index. Initialize with -1, such that empty bins, and the merged bins pointing to a removed IBF will have a pointer to -1.
    int counter = 0;// Initialize the result vector
    for (size_t i = 0; i < index.ibf().ibf_vector.size(); i++) {
        if (std::find(indices_to_remove.begin(), indices_to_remove.end(), i) == indices_to_remove.end()) {  // If the current element is not in indices_to_remove.
            indices_map[i] = counter; // Add it to the result vector, such that indices_map[i] = counter. Do not use .push_back, because that will make the list not long enough
            counter += 1;
        }
    };

    // Remove vectors of indices of subindex datastructures like next_ibf and previous_ibf
    remove_indices(indices_to_remove, index.ibf().ibf_vector);
    remove_indices(indices_to_remove, index.ibf().next_ibf_id);
    remove_indices(indices_to_remove, index.ibf().previous_ibf_id);
    remove_indices(indices_to_remove, index.ibf().fpr_table);
    remove_indices(indices_to_remove, index.ibf().occupancy_table);
    remove_indices(indices_to_remove, index.ibf().user_bins.ibf_bin_to_filename_position);
    for (size_t ibf_idx{0}; ibf_idx < index.ibf().next_ibf_id.size(); ibf_idx++) {     // Replace the (remaining) indices that have to be replaced.
        for (size_t i{0}; i < index.ibf().next_ibf_id[ibf_idx].size(); ++i){
            auto & next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][i];
            if (next_ibf_idx >=0){
                if (static_cast<size_t>(next_ibf_idx) >= indices_map.size())
                    next_ibf_idx = -1;
                else
                //assert(static_cast<size_t>(next_ibf_idx) < indices_map.size()); // this assertion prevents an out of range segmentation fault
                    next_ibf_idx = indices_map.at(next_ibf_idx); // TODO check another time if this goes as should
                // easiest but not the most neat solution is to convert it to a -1.
                //                 assert(static_cast<size_t>(next_ibf_idx) < indices_map.size()); // this assertion prevents an out of range segmentation fault
                //                next_ibf_idx = indices_map[next_ibf_idx];
            }
        }
    }
    for (size_t ibf_idx = 0; ibf_idx < index.ibf().previous_ibf_id.size(); ++ibf_idx){
        auto & previous_ibf_idx = std::get<0>(index.ibf().previous_ibf_id[ibf_idx]);
        previous_ibf_idx = indices_map[previous_ibf_idx];
    }
}

/*!\brief Merges the original HIBF with the pruned subtree, with the rebuild subtree.
 * \details One should remove the IBFs in the original index which were part of the subtree that had to be rebuild.
 * When doing some sort of splitting on the merged bins, run this function twice, once for each subindex.
 * \param[in|out] index the original HIBF
 * \param[in] subindex the HIBF subtree that has been rebuild
 * \param[in] ibf_idx the index of the IBF where the subtree needs to be removed. (irrelevant?)
 * \param[in] index_tuple the index of the IBF and the bin index of the MB to which the subtree needs to be attached.
 * \author Myrthe Willemsen
 */
void attach_subindex(raptor_index<index_structure::hibf> & index,
                     raptor_index<index_structure::hibf> & subindex,
                     std::tuple<size_t, size_t> index_tuple){
    // Add new rows representing the subindex.
    size_t ibf_count_before_appending = index.ibf().ibf_count();
    for (size_t ibf_idx{0}; ibf_idx < subindex.ibf().next_ibf_id.size(); ++ibf_idx){ // Add the size of the `index`, in number of IBFs, to all IBF indices in subindex's next_ibf_id.
        for (size_t bin_idx{0}; bin_idx < subindex.ibf().next_ibf_id[ibf_idx].size(); ++bin_idx){
             subindex.ibf().next_ibf_id[ibf_idx][bin_idx] += ibf_count_before_appending;
        }
        std::get<0>(subindex.ibf().previous_ibf_id[ibf_idx]) += ibf_count_before_appending; // Add the size of the `index`, in number of IBFs, to the IBF indices present in previous_ibf_id of the subindex.
    }

    auto append_to_vector = [] (auto & index_vector, auto & subindex_vector){
        index_vector.insert(index_vector.end(), subindex_vector.begin(), subindex_vector.end());
    };
    append_to_vector(index.ibf().ibf_vector, subindex.ibf().ibf_vector); // append the vector of the new index to the original index.
    append_to_vector(index.ibf().next_ibf_id, subindex.ibf().next_ibf_id);
    append_to_vector(index.ibf().previous_ibf_id, subindex.ibf().previous_ibf_id);
    append_to_vector(index.ibf().fpr_table, subindex.ibf().fpr_table);
    append_to_vector(index.ibf().occupancy_table, subindex.ibf().occupancy_table);

    // Map the filename indices of the subindex to those of the original index. Then you do not additionally need to append the filenames vectors themselves.
    size_t subindex_n_files = subindex.ibf().user_bins.user_bin_filenames.size();
    std::vector<int> indices_map; indices_map.resize(subindex_n_files); // The map has the size in number of files in the subindex.
    std::ranges::fill(indices_map, -1); // initialize with -1, such that empty bins, and the merged bins pointing to a removed IBF will have a pointer to -1.
    for (size_t filename_idx{0}; filename_idx < subindex_n_files; filename_idx++){
        auto filename = subindex.ibf().user_bins.user_bin_filenames[filename_idx];
        if (std::filesystem::path filename_as_path{filename};
        filename_as_path.extension() != ".empty_bin" and filename != ""){
            subindex.ibf().user_bins.user_bin_filenames[filename_idx];
            indices_map[filename_idx] = (index.ibf().user_bins.filename_to_idx.at(filename));
        }
    } // mapping works as follows: filename_idx --> idx_to_filename of subindex --> filename --> filename_to_idx of original index --> idx

    auto & bin_to_file = subindex.ibf().user_bins.ibf_bin_to_filename_position;
    for (size_t ibf_idx{0}; ibf_idx < bin_to_file.size(); ++ibf_idx){ // Add the size of the `index`, in number of IBFs, to all IBF indices in subindex's next_ibf_id.
        for (size_t bin_idx{0}; bin_idx < bin_to_file[ibf_idx].size(); ++bin_idx){
            if (bin_to_file[ibf_idx][bin_idx] != -1){ // if it is not an empty bin .
                assert(bin_to_file[ibf_idx][bin_idx] < indices_map.size());
                bin_to_file[ibf_idx][bin_idx] = indices_map[bin_to_file[ibf_idx][bin_idx]];
            }
        }
    }
    append_to_vector(index.ibf().user_bins.ibf_bin_to_filename_position, bin_to_file);

    // Update the indices in one entry of the supporting tables, where are subindex must be attached, such that they refer to our new subindex.
    auto ibf_idx = std::get<0>(index_tuple);
    auto bin_idx = std::get<1>(index_tuple);
    index.ibf().next_ibf_id[ibf_idx][bin_idx] = ibf_count_before_appending;
    index.ibf().previous_ibf_id[ibf_count_before_appending] = std::make_tuple(ibf_idx, bin_idx);
    index.ibf().user_bins.ibf_bin_to_filename_position[ibf_idx][bin_idx] = -1; // filename index is set to -1 (filename is undefined) for merged bins.
}

} // namespace raptor

