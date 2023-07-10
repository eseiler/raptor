#include <raptor/update/insertions.hpp>
#include <raptor/update/rebuild.hpp>
#include <raptor/build/hibf/compute_kmers.hpp>
#include <raptor/build/hibf/insert_into_ibf.hpp>
#include <numeric>
#include <chopper/sketch/hyperloglog.hpp>
#include "chopper/sketch/output.hpp"
#include "chopper/sketch/toolbox.hpp"
#include "chopper/sketch/execute.hpp"
#include "chopper/next_multiple_of_64.hpp"

namespace raptor
{

size_t new_bin_count(size_t number_of_bins, double eb_fraction, size_t ibf_bin_count){
    return chopper::next_multiple_of_64(std::max((size_t) std::round((1+eb_fraction)
        *ibf_bin_count), ibf_bin_count + number_of_bins));
};

/*!\brief Finds a location in the HIBF for a new UB.
 * \details The algorithm finds a suitable location in the existing HIBF for a new UB, by .
 * 1) looking for a proper IBF with one of the three algorithms
 * 2) calculating the number of TBs needed to insert the UB in this IBF
 * 3) finding the bin_idx in the
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \return ibf_idx, bin_idx, number_of_bins An index triple of the index of IBF, start index of the technical bins, and the number of bins
 * \author Myrthe Willemsen
 */

std::tuple <uint64_t, uint64_t, uint16_t> get_location(size_t kmer_count, robin_hood::unordered_flat_set<size_t> &kmers,
                                                       raptor_index<index_structure::hibf> & index,
                                                       update_arguments const & update_arguments){
    size_t root_idx = 0; size_t ibf_idx=-1;
    if (update_arguments.ibf_selection_method == "find_ibf_idx_traverse_by_similarity")
        ibf_idx = find_ibf_idx_traverse_by_similarity(kmers, index);
    else if (update_arguments.ibf_selection_method == "find_ibf_idx_ibf_size")
        ibf_idx = find_ibf_idx_ibf_size(kmers.size(), index);
    else if (update_arguments.ibf_selection_method == "find_ibf_size_splitting")
        ibf_idx = find_ibf_size_splitting(kmers.size(), index, update_arguments);
    else if (update_arguments.ibf_selection_method == "find_ibf_idx_traverse_by_fpr"){
        size_t ibf_idx_by_size = find_ibf_idx_traverse_by_fpr(kmer_count, index, root_idx);
        if (update_arguments.tmax_condition){
            ibf_idx = find_ibf_idx_traverse_by_fpr_tmax(kmer_count, index, update_arguments, ibf_idx_by_size);
            if (ibf_idx == -1)
                ibf_idx =  ibf_idx_by_size; // in fact, we should do here a partial rebuild, adding the new filename to the subtree. Currently a rebuild is triggered because we resize beyond tmax.  TODO
        }
    }
    else
        std::cout << "No correct IBF selection method was entered";
    size_t number_of_bins = 1; // calculate number of user bins needed.
    if (index.ibf().ibf_max_kmers(ibf_idx) < (size_t) kmer_count){ // when using the find_ibf_idx_ibf_size, it segfaults in the bin_size function, when size_t bin_size = ibf.bin_size() in ibf_max_kmers
        // Only if we are at the root we might have to split bins. Delete ibf_idx==0, if using similarity method
        number_of_bins = index.ibf().number_of_bins(ibf_idx, (int) kmer_count);       // calculate among how many bins we should split
    }
    uint64_t bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, number_of_bins);
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    if (bin_idx == ibf_bin_count){ // current solution TODO
        size_t new_ibf_bin_count = new_bin_count(number_of_bins, update_arguments.empty_bin_percentage, ibf_bin_count);
        std::cout << "Resize the IBF at index: " << ibf_idx << "\n" << std::flush;
        index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
    }

    return {ibf_idx, bin_idx, number_of_bins};
    }


//!\brief  Compute and store sketches for the new UBs using the choppper library
void update_sketch(std::vector<std::string> filename, update_arguments update_arguments,
                raptor_index<index_structure::hibf> & index, bool existing_sketch=false){

    if (update_arguments.sketch_directory != ""){ // if sketches are used, then update the sketch of this UB.
        std::vector<chopper::sketch::hyperloglog> sketches;
        if (existing_sketch){
            chopper::sketch::toolbox::read_hll_files_into(update_arguments.sketch_directory, filename, sketches); // instead of hll_dir, use update_arguments.sketch_directory
        }
        chopper::configuration layout_arguments = layout_config(index, update_arguments); // Create the arguments to run the layout algorithm with.
        chopper::sketch::execute(layout_arguments, filename, sketches); // make sure that config.precomputed_files = false.
    }
}
/*!\brief Inserts a UB in the assigned TBs and its parent MBs in higher level IBFs
 * \details The algorithm inserts the UB in the leave bins (LBs) and parent merged bins
 * iteratively until reaching the root ibf at the top level. For each level the function
 * `insert_into_ibf` will be called. This function updates the `rebuild_index_tuple`, indicating
 * whether part of the index needs to be rebuild. This is initialized with the ibf_idx being the HIBF size.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index_triple (ibf_idx, bin_idx, number_of_bins) indicating where the UB needs to be inserted. Not by reference, because the index_triple is modified.
 * \param[in] index the original HIBF
 * \return rebuild_index_tuple containing (ibf_idx, bin_idx). If the ibf_idx is lower than the total number of IBFs,
 * then the respective IBF needs to be rebuild.
 * \author Myrthe Willemsen
 */
std::tuple <uint64_t, uint64_t> insert_tb_and_parents(robin_hood::unordered_flat_set<size_t> & kmers,
                                                      std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                                                      raptor_index<index_structure::hibf> & index){
    std::tuple <uint64_t, uint64_t> rebuild_index_tuple = std::make_tuple(index.ibf().ibf_vector.size(), 0);
    while (true){
        insert_into_ibf(kmers, index_triple, index, rebuild_index_tuple);
        if (std::get<0>(index_triple) == 0) break;
        auto index_tuple = index.ibf().previous_ibf_id[std::get<0>(index_triple)]; // update index triple:
        index_triple = std::make_tuple(std::get<0>(index_tuple), std::get<1>(index_tuple), 1); // number of bins will be 1 for the merged bins (i assume)
        }
    return rebuild_index_tuple;
    }

/*!\brief Insert (multiple) new UBs
 * \details The algorithm interts new UBs. At last it computes and stores a sketch for the new UB and stores this.
 * The user is supposed to add the filename to the file containing the bin paths (all_bins_path) himself, if desired.
 * \param[in] update_arguments
 * \param[in] index the original HIBF
 * \author Myrthe Willemsen
 */
void insert_ubs(update_arguments const & update_arguments,
                raptor_index<index_structure::hibf> & index){
    for  (auto &filename: update_arguments.bin_path){ // Loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        std::cout << "Exists filename method \n";
        if (index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin, inserts it if it does not exist yet.
            std::cout << "The user bin ... that you want to insert does already exist. Please use the --insert-sequences option";
        }else{
            std::string filename_cast = filename[0];
            std::cout << "Inserting bin: " << filename_cast << "\n";
            robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize kmers.
            assert(update_arguments.shape.size() > 0);
            raptor::hibf::compute_kmers(kmers, update_arguments, filename);
            size_t kmer_count = kmers.size();
            std::tuple <uint64_t, uint64_t, uint16_t> index_triple = get_location(kmer_count, kmers, index, update_arguments);  //  index_triple; bin_idx, ibf_idx, number_of_bins
            std::tuple <uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers,  index_triple, index);

            index.ibf().user_bins.update_filename_indices(filename_cast, index_triple);
            update_sketch(filename, update_arguments, index, false);             //Compute and store sketches for the new UBs

            if (std::get<0>(rebuild_index_tuple) < index.ibf().ibf_vector.size()) // If the ibf_idx is lower than the total number of IBFs (the initilization value), then the respective IBF needs to be rebuild.
                partial_rebuild(rebuild_index_tuple, index, update_arguments);
            size_t const ibf_idx = std::get<0>(index_triple);

            check_tmax_rebuild(index, ibf_idx, update_arguments); // check if rebuilding is needed because of exceeding the tmax
        }
    }
}


/*!\brief splits UB content over more TBs
 * \details The algorithm determines the new number of TBs to be split over, using the percentage of empty bins.
 * \param[in] index_triple
 * \param[in] kmers the set of new kmers (the existing k-mers are added to these)
 * \param[in] update_arguments
 * \param[in] index the original HIBF
 * \author Myrthe Willemsen
 */
void split_user_bin(std::tuple <uint64_t, uint64_t, uint16_t> index_triple,
                    robin_hood::unordered_flat_set<size_t> & kmers,
                    update_arguments const & update_arguments,
                    raptor_index<index_structure::hibf> & index,
                    std::string filename){
    size_t ibf_idx = std::get<0>(index_triple);
    size_t start_bin_idx = std::get<1>(index_triple);
    size_t number_of_bins = std::get<2>(index_triple);
    index.ibf().delete_tbs(ibf_idx, start_bin_idx, number_of_bins);
    raptor::hibf::compute_kmers(kmers, update_arguments, std::vector{filename}); // Get existing k-mers. Those are appended to the k-mer array that is already loaded.
    size_t new_number_of_bins = std::ceil((int) number_of_bins * (1 + update_arguments.empty_bin_percentage));  // new number of bins * empty_bin_percentage
    new_number_of_bins = std::max(new_number_of_bins, index.ibf().number_of_bins(ibf_idx, (int) kmers.size()));       // make sure that the new number of bins is also sufficient to store the additional k-mers.
    size_t new_start_bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, new_number_of_bins); // Find empty bins.
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    if (new_start_bin_idx == ibf_bin_count){ // current solution. Alternatively, (1) calculate the new bin count of the IBF (2) check if this is larger than the tmax (3 - no) insert (3 yes) insert only in parents en rebuild subtree
        size_t new_ibf_bin_count = new_bin_count(new_number_of_bins, update_arguments.empty_bin_percentage, ibf_bin_count);
        std::cout << "Resize the IBF at index: " << ibf_idx << "\n" << std::flush;
        index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
    }
    std::tuple <uint64_t, uint64_t, uint16_t> new_index_triple = std::make_tuple(ibf_idx, new_start_bin_idx, new_number_of_bins);
    insert_tb_and_parents(kmers, new_index_triple, index); // Insert the k-mer content of the user bin again.
    index.ibf().user_bins.update_filename_indices(filename, new_index_triple); // update additional datastructures.
    check_tmax_rebuild(index, ibf_idx, update_arguments); // check if rebuilding is needed because of exceeding the tmax
}

/*!\brief Inserts sequences in existing UBs
 * \details The algorithm inserts new sequence content in the technical bins of an existing sample,
 * as well as its parent merged bins.
 * If sketches are used, then also the sketch of this UB is updated.
 * \guideline The user can provide a file specifically containing the part of the sequence that should be added, not the whole sequence.
 * By default this file ends with "_insertsequences", but can be set to a different appendix.
 * I.e. if you want to insert to the sample "bin_00.fasta", then the file containing the new sequences should be called "bin_00_insertsequences.fasta"
 * If this is not possible, the user can provide the
 * I should update hyperloglog sketches, wheres the user: updates sequence files yourself., e.g. using cat.
 * \param[in] kmers the set of kmers to be stored
 * \param[in] index the original HIBF
 * \author Myrthe Willemsen
 */
void insert_sequences(update_arguments const & update_arguments, raptor_index<index_structure::hibf> & index){
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        if (not index.ibf().user_bins.exists_filename(filename[0])){ // Find location of existing user bin. Note that this function inserts it if it does not exist yet.
            std::cout << "The user bin ... that you want to insert to does not exist. If you want to add a new user bin, use the flag -insert-UB";
        }else{
                robin_hood::unordered_flat_set<size_t> kmers{}; // Initialize k-mers.
                std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
                assert(std::get<2>(index_triple)); // assert that the num_user_bins is not 0.
                std::string filename_new_sequences = filename[0].substr(0, filename[0].find_last_of('.')) + update_arguments.insert_sequence_appendix + filename[0].substr(filename[0].find_last_of('.'));
                raptor::hibf::compute_kmers(kmers, update_arguments, std::vector{filename_new_sequences});
                std::tuple <uint64_t, uint64_t> rebuild_index_tuple = insert_tb_and_parents(kmers, index_triple, index);
                update_sketch(filename, update_arguments, index, true);
                if (std::get<0>(rebuild_index_tuple) < index.ibf().ibf_vector.size()){ // If the ibf_idx is lower than the total number of IBFs (the initilization value), then the respective IBF needs to be rebuild.
                    if (std::get<0>(rebuild_index_tuple) != std::get<0>(index_triple))
                        partial_rebuild(rebuild_index_tuple, index, update_arguments, 2);
                    else split_user_bin(index_triple, kmers, update_arguments, index, filename[0]); // If the user bin reaches the FPRmax, then it can simply be split over more technical bins.
                }

        }
    }
}



/*!\brief Delete a multiple user bins from the index.
 * \param[in] update_arguments
 * \param[in] index The HIBF.
 * \guide deleting the actual sequence file and from all_bins_path is to the user.
 * \author Myrthe Willemsen
 */
void delete_ubs(update_arguments const & update_arguments,
                  raptor_index<index_structure::hibf> & index){
    for  (auto &filename: update_arguments.bin_path){ // loop over new bins, using arguments.bin_path, as created in parse_bin_path(arguments) in upgrade_parsing.cpp
        delete_ub(filename, index); // delete a single user bin from the index.

        if (update_arguments.sketch_directory != ""){ // if sketches are used, then delete the sketch of this UB.
            std::filesystem::path path = update_arguments.sketch_directory / std::filesystem::path(filename[0]).stem();
            path += ".hll";
            std::filesystem::remove(path);
        }
    }
}

/*!\brief Delete a single user bin from the index.
 * \param[in] filename filename of the user bin to be removed.
 * \param[in] index The HIBF.
 * \details Delete a single UB by deleting one ore more leaf bins,
 * with the help of the the `delete_tbs` function.
 * \author Myrthe Willemsen
 */
void delete_ub(std::vector<std::string> const & filename,
                    raptor_index<index_structure::hibf> & index){

    if (not index.ibf().user_bins.exists_filename(filename[0])) // first find index
    {
        std::cout << "Warning: the user bin you want to delete is not present in the HIBF: "; // + filename; //--> if not: return error , make sure by doing this, we dont accidently add the filename..
    }else{
        std::tuple <uint64_t, uint64_t, uint16_t> index_triple = index.ibf().user_bins.find_filename(filename[0]);
        size_t const ibf_idx = std::get<0>(index_triple); // Create an empty UB
        size_t const start_bin_idx = std::get<1>(index_triple);
        size_t const number_of_bins = std::get<2>(index_triple);
        index.ibf().delete_tbs(ibf_idx, start_bin_idx, number_of_bins);
        index.ibf().user_bins.delete_filename(filename[0]);  // Update filename tables. Even if the UB did not exist, it might have been added through the STL's .find() function.
    }
}





/*!\brief Finds an appropiate IBF for a UB insertion based on FPR.
 * \param[in] kmer_count The number of k-mers to be inserted.
 * \param[in] ibf_idx The IBF index of the current IBF. Should be the root (0) at the start of the traversal.
 * \param[in] index The HIBF.
 * \return the IBF index which is best appropiate based on the FPR, using the FPR table.
 * \details the merged bin with the lowest false positive rate is selected.
 * The algorithm recurses until an IBF's bin size is too small to accommodate the new UB.
 * \author Myrthe Willemsen
 */
size_t find_ibf_idx_traverse_by_fpr(size_t & kmer_count, raptor_index<index_structure::hibf> & index, size_t ibf_idx =0){ //default is root_idx =0, where we start the search.
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF
    if (index.ibf().ibf_max_kmers(ibf_idx) > kmer_count){ // kmer-capacity of IBF > bin size new UB, go down if possible. Instead of maximal capcity, you can calculate the optimal kmer_ size.
        size_t best_mb_idx = ibf.bin_count(); double best_fpr = 1; // initialize the best idx outside of the ibf, such that we can use this after the loop to check if a MB was found.
         for (size_t bin_idx=0; bin_idx < ibf.bin_count(); ++bin_idx){ //loop over bins to find the bext merged bin
            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
                auto fpr = index.ibf().get_fpr(ibf_idx, bin_idx);
                if (fpr < best_fpr){
                    best_fpr = fpr;
                    best_mb_idx = bin_idx;
                }
            }
         }
         if (best_mb_idx >= ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
             return ibf_idx;
         }else{
             auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //next_ibf_id[ibf_id_high].size()
             return (find_ibf_idx_traverse_by_fpr(kmer_count, index, next_ibf_idx));
         }
    }else{ // kmer-capacity of IBF < bin size new UB, go up if possible
        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
    }
    }

/*!\brief Traverse the tree further down until an IBF is found with sufficient empty bins.
 * \param[in] ibf_idx of the subtree
 * \param[in] number_of_bins the number of bins needed to store the UB in question in this particular IBF.
 * \param[in] index The HIBF.
 * \return The IBF in the subtree that has sufficient empty bins
 * \details If the the `tmax` is used as a threshold on the number of technical bins, then a good insertion method should
 * take this into account.
 * At each IBF we will use the function number_of_bins to calculate the required number of bins, and then the
 * find_empty_bin_idx function to check if there are sufficient empty bins.
 * \author Myrthe Willemsen
 */
 size_t find_ibf_idx_traverse_by_fpr_tmax(size_t & kmer_count, raptor_index<index_structure::hibf> & index,
                                          update_arguments const & update_arguments,size_t ibf_idx){
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF
    auto number_of_bins = index.ibf().number_of_bins(ibf_idx, (int) kmer_count);       // calculate among how many bins we should split
    size_t start_bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, number_of_bins); // Find empty bins.

    if (start_bin_idx != ibf.bin_count())
        return(ibf_idx); // if we do find empty bins, just return this IBF.
    else if (start_bin_idx == ibf.bin_count()){ // If we can find no empty bins in the IBF, traverse down.
        size_t best_mb_idx = ibf.bin_count(); double best_fpr = 1; // initialize the best idx outside of the ibf, such that we can use this after the loop to check if a MB was found.
         for (size_t bin_idx=0; bin_idx < ibf.bin_count(); ++bin_idx){ //loop over bins to find the bext merged bin
            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
                auto fpr = index.ibf().get_fpr(ibf_idx, bin_idx);
                if (fpr < best_fpr){
                    best_fpr = fpr;
                    best_mb_idx = bin_idx;
                }
            }
         }
         if (best_mb_idx >= ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
             return -1; //if no IBFs have been found in this subtree, then insert somewhere up the subtree, and rebuild? ; for now, insert in the original level, which triggers a rebuild itself. This can be improved by just inserting the sequence content in the merged bins up the tree, and then rebuilding without inserting the file beforehand.
         }else{
             auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //enter the next IBF.
             return (find_ibf_idx_traverse_by_fpr_tmax(kmer_count, index, update_arguments, next_ibf_idx));
         }
    }
}


/*!\brief Finds empty bins (EBs) within a certain IBF, where the new UB can be inserted.
 * \param[in] ibf_idx the IBF in which EBs need to be found.
 * \param[in] number_of_bins the number of bins needed to store the UB in question in this particular IBF.
 * \param[in] index The HIBF.
 * \return The starting index of the empty TBs where the new bin can be inserted.
 * \details Try inserting in the first EB encountered. Check if there are sufficient adjacent empty bins.
 * Note: this could be improved using empty bin data structure and rank operation. Also jump table search would help
 * to speed up this method.
 * If there is no empty bin, the IBF must be resized.
 * \author Myrthe Willemsen
 */
uint64_t find_empty_bin_idx(raptor_index<index_structure::hibf> & index, size_t ibf_idx,
                            update_arguments const & update_arguments, size_t number_of_bins)
    {
    size_t ibf_bin_count = index.ibf().ibf_vector[ibf_idx].bin_count();
    size_t bin_idx{0}; // The variable is initialized outside the for loop, such that afterwards it can still be used.
    for (; bin_idx + number_of_bins < ibf_bin_count; bin_idx++){ // This could be implemented more efficiently.
        if (std::reduce(&index.ibf().occupancy_table[ibf_idx][bin_idx], // reduce sums over all values in the range
                        &index.ibf().occupancy_table[ibf_idx][bin_idx + number_of_bins])==0 //this range is empty, so a good location has been found. Std reduce should go from [bin_idx] to [bin_idx_end - 1]
                        or std::reduce(&index.ibf().occupancy_table[ibf_idx][bin_idx], // reduce sums over all values in the range
                        &index.ibf().occupancy_table[ibf_idx][bin_idx + number_of_bins])==number_of_bins){ // sometimes empty bins are filled with 1, a problem not resolved. The sum is then number_of_bins.
            if (bin_idx + number_of_bins - 1 < ibf_bin_count) return bin_idx; // bin_idx + number_of_bins - 1 is the last index that would be occupied.

    } // If nothing has been returned, no appropriate empty bin has been found and the bin idx will be the size of the IBF,
    bin_idx = ibf_bin_count; // then the IBF must be resized.
    size_t new_ibf_bin_count = chopper::next_multiple_of_64(std::max((size_t) std::round((1+update_arguments.empty_bin_percentage)
            *ibf_bin_count), ibf_bin_count + number_of_bins));
    assert(new_ibf_bin_count > ibf_bin_count); // make sure that the new bin count is larger than the current IBF size.
    if (new_ibf_bin_count <= index.ibf().t_max){
        std::cout << "Resize the IBF at index: " << ibf_idx << "\n" << std::flush;
        index.ibf().resize_ibf(ibf_idx, new_ibf_bin_count);
        assert(index.ibf().occupancy_table[ibf_idx][bin_idx] == 0); // check if the newly added bins are really 0 .
    }
    return bin_idx;
}

/*!\brief Finds an appropiate IBF for a UB insertion based on the number of k-mers and IBF sizes.
 * \param[in] kmer_count The number of k-mers to be inserted.
 * \param[in] index The HIBF.
 * \return the IBF index to insert the new UB in.
 * \details A second algorithm picks the IBF that has the smallest size able to store the UB without it being split.
 * It finds this IBF by doing a binary search on a sorted array with IBF bin sizes in terms of the number of
 * $k$-mers that they can maximally store, and the corresponding IBF indexes. The search is logarithmic in the time of the number of IBFs,
 * \author Myrthe Willemsen
 */
size_t find_ibf_idx_ibf_size(size_t kmer_count, raptor_index<index_structure::hibf> & index){
    auto & array = index.ibf().ibf_sizes;
        int low = 0;
        int high = array.size()-1;
        while (low <= high) {
            int mid = (low + high) >> 1;
            assert(mid < array.size());
            if (std::get<0>(array[mid]) < kmer_count)
                {low = mid + 1;}
            else if (std::get<0>(array[mid]) > kmer_count)
                {high = mid - 1;}
            else if (std::get<0>(array[mid]) == kmer_count)
                {return std::get<1>(array[mid]);} // exact kmer_count found
        }
        low = std::min(low, static_cast<int>(array.size())-1); // low = mid + 1, so it may happen that low equals the array.size.
        assert(low < array.size());
        return std::get<1>(array[low]);
    }


/*!\brief Finds an appropiate IBF for a UB insertion based on the number of k-mers and IBF sizes.
 * \param[in] kmer_count The number of k-mers to be inserted.
 * \param[in] index The HIBF.
 * \return the IBF index to insert the new UB in.
 * \details
 * \author Myrthe Willemsen
 */
size_t find_ibf_size_splitting(size_t kmer_count, raptor_index<index_structure::hibf> & index,
                               update_arguments const update_arguments){
    auto & array = index.ibf().ibf_sizes;
        int low = 0;
        int high = array.size()-1;
        while (low <= high) {
            int mid = (low + high) >> 1;
            assert(mid < array.size());
            if (std::get<0>(array[mid]) < kmer_count)
                {low = mid + 1;}
            else if (std::get<0>(array[mid]) > kmer_count)
                {high = mid - 1;}
            else if (std::get<0>(array[mid]) == kmer_count)
                {low = mid; break;} // exact kmer_count found
        }
        // SPLITTING BELOW
        // If the found IBF ha sno empty bin, then check for each of the IBFs with a smaller size whether they have
        // sufficient empty bins to insert the new user bin. If so, return the IBF index. If not, continue with the next section.
        low = std::min(low, static_cast<int>(array.size())-1); // low = mid + 1, so it may happen that low equals the array.size.
        assert(low < array.size());
        auto low_perfect = low;
        while (low >= 0){
            auto ibf_idx = std::get<1>(array[low]);
            auto& ibf = index.ibf().ibf_vector[ibf_idx]; //  select the IBF
            auto number_of_bins = index.ibf().number_of_bins(ibf_idx, (int) kmer_count);       // calculate among how many bins we should split
            size_t start_bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, number_of_bins); // Find empty bins.
            if (start_bin_idx + number_of_bins <= ibf.bin_count())
                return(ibf_idx); // if we do find empty bins, just return this IBF.
            else low -= 1;
        }

        // INSERTING IN LARGER IBFs
        // Check for the parent ibf above if there is still an empty bin such that a partial rebuild can be done,
        // If it has no space for empty bins before reaching the tmax, then check the IBFs with sizes between the original IBF and the parent IBF.
        // Should any of them have an empty bin, then insert the user bin there.
        // If not, then repeat the process for the parent of the parent.
        // This way, empty bins will be filled up before the next full rebuild.
        low = low_perfect;
        size_t ibf_idx = std::get<1>(array[low]);

        while (ibf_idx){ // while the ibf is not the root.
            size_t parent_ibf_idx = std::get<0>(index.ibf().previous_ibf_id[ibf_idx]);
            auto ibf_parent = index.ibf().ibf_vector[parent_ibf_idx];
            if (static_cast<size_t>(find_empty_bin_idx(index, parent_ibf_idx, update_arguments, 1))
                != ibf_parent.bin_count()){ // the parent has at least one empty bin, which should accomodate a partial rebuild.
                    return ibf_idx; // this should trigger a resize and partial rebuild down the stream.
            }else{
                auto size_parent = ibf_parent.bin_size();
                //low = low_perfect;
                while (low < array.size() and std::get<1>(array[low]) < size_parent ){ // and thereby low will be lower than the maximum value of the array
                    auto ibf_idx = std::get<1>(array[low]);
                    size_t start_bin_idx = find_empty_bin_idx(index, ibf_idx, update_arguments, 1); // Find empty bins.
                    if (start_bin_idx != index.ibf().ibf_vector[ibf_idx].bin_count())
                        return(ibf_idx); // if we do find empty bins, just return this IBF.
                    else
                        low += 1;
                }
            }
            size_t ibf_idx = parent_ibf_idx; // would be equal to using  ibf_idx= std::get<1>(array[low]);
        }
        return ibf_idx; // a full rebuild will be triggered.
    }



/*!\brief Finds an appropiate IBF for a UB insertion based on sequence similarity
 * \param[in] kmers The kmers to be inserted.
 * \param[in] ibf_idx The IBF index of the current IBF. Should be the root (0) at the start of the traversal.
 * \param[in] index The HIBF.
 * \return ibf_idx of the IBF to insert the new UB in.
 * \details the merged bin with the maximal similarity is selected.
 * The algorithm recurses until an IBF's bin size is too small to accommodate the new UB.
 * The similarities to all TBs in each traversed IBF are approximated by querying a percentage of the k-mers from the new UBs.
 * This can be efficiently done using the query me
 * thod that is already in place.
 * \author Myrthe Willemsen
 */
size_t find_ibf_idx_traverse_by_similarity(robin_hood::unordered_flat_set<size_t> & kmers, raptor_index<index_structure::hibf> & index, size_t ibf_idx){ //default is root_idx =0, where we start the search.
    auto& ibf = index.ibf().ibf_vector[ibf_idx]; //select the IBF
    if (index.ibf().ibf_max_kmers(ibf_idx) > kmers.size()){ // kmer-capacity of IBF > bin size new UB, go down if possible. Instead of maximal capcity, you can calculate the optimal kmer_ size.
        auto agent = ibf.template counting_agent<uint16_t>();
        int num_kmers_to_sample = 0.01*kmers.size(); // Create a new unordered_flat_set to store the sampled kmers
        std::vector<size_t> sampled_kmers;
        std::ranges::sample(kmers, std::back_inserter(sampled_kmers), num_kmers_to_sample, std::mt19937 {std::random_device {}()});

        auto & result = agent.bulk_count(sampled_kmers); // count occurrences of the kmers in each of the bins in the current IBF.
        size_t best_mb_idx = ibf.bin_count(); int best_similarity = -1; // initialize the best index outside of the ibf, such that we can use this after the loop to check if a MB was found.
        for (size_t bin_idx=0; bin_idx < ibf.bin_count(); ++bin_idx){ // loop over bins to find the next merged bin
            if (index.ibf().is_merged_bin(ibf_idx, bin_idx)){
                auto similarity = result[bin_idx];
                if (similarity > best_similarity){
                    best_similarity = similarity;
                    best_mb_idx = bin_idx;
                }
            }
        }
        if (best_mb_idx >= ibf.bin_count()){ //no merged bins, only leaf bins exist on this level.
            return ibf_idx;
        }else{
            auto next_ibf_idx = index.ibf().next_ibf_id[ibf_idx][best_mb_idx]; //next_ibf_id[ibf_id_high].size()
            return (find_ibf_idx_traverse_by_similarity(kmers, index, next_ibf_idx));
        }
    }else{ // if the kmer-capacity of the IBF is lower than the required bin size for the new UB, then go up if possible
        return(std::get<0>(index.ibf().previous_ibf_id[ibf_idx])); //ibf idx of merged bin a level up. If it is a root, it will automatically return the root index = 0
    }
}








} // end namespace
