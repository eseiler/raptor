#include <raptor/search/search_single.hpp>
#include <raptor/update/load_hibf.hpp>

#include <raptor/build/store_index.hpp>

//


namespace raptor
{

/*!\brief Loads the stored HIBF into memory.
 * \remark does not yet work for compressed indexes.
 * \param[out] index the HIBF
 * \param[in] update_arguments configuration object with parameters required for calling an update operation
 * \author Myrthe Willemsen
 */
void load_hibf(raptor_index<index_structure::hibf> & index, update_arguments & arguments) // perhaps better to have index as in and output of this function, because now it is calling the update function within.
{   if (not arguments.compressed){
    double index_io_time{0.0};
    load_index(index, arguments, index_io_time);
    }else{ // try decompressing a IBF
        auto some_compressed_ibf = index.ibf().ibf_vector[0];
        seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf_uncompressed{some_compressed_ibf};
    }
    arguments.shape = index.shape();
    arguments.window_size = index.window_size();
    assert(arguments.shape.size() > 0);
    assert(arguments.shape.count() > 0);
    return;
}

} // end namespace

