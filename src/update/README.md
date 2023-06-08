##Updating the HIBF

This readme file provides instructions on updating the Dynamic Hierarchical Interleaved Bloom Filter (HIBF), grouped by different update operations.

In general, the fraction of empty bins (`empty_bin_percentage`) determines the fraction of resizing during updates. 
The default value is set to `0.5`, but it can be adjusted as needed.


### User Bin Modifications
The user can update the content of a user bin by the flag `--insert-sequences`. The user should provide a file containing only the specific part of the sequence that needs to be added. The filename containing the new sequence content is the original filename, ending on
a specific appendix, which can be input by `--insert_sequence_appendix`. By default the ending is `"_insertsequences"`. It is the responsibility of the user to update the fasta files themselves, 
such as by using the concatenate `cat` command. The HyperLogLog sketches are automatically updated during the index update process. 
Note that sequence deletions are not supported. 

### User Bin Insertions 
If we want to insert specific samples from the index, we provide a list of their filenames through the paremeter `--bins` and use the flag `--insert-UBs`.
We can specify the method for finding the location for the new user bin, 
by the parameter `--ibf_selection_method`, with the options `find_ibf_idx_traverse_by_similarity` (by default), `find_ibf_idx_ibf_size` or `find_ibf_idx_traverse_by_fpr`.


### User Bin Deletions
If we want to delete specific samples from the index, we provide a list of their filenames through the paremeter `--bins` and use the flag `--delete-UBs`. 
The content of the corresponding technical bins IBF (Interleaved Bloom Filter) vector will be cleared accordingly.
Please note that if a significant portion of the bins is removed, it might be advisable to rebuild the entire data structure. 
Unlike a dynamic array, the HIBF does not automatically shrink when content disappears, 
as it is not expected to occur frequently. In that case one should ensure that the file containing the bin paths (`bins`) used for computing the layout is updated by the user when performing updates or rebuilding the index.
Deleting Original FASTA Files: If desired, the user should manually delete the original FASTA file of the deleted samples, as otherwise it will keep occupying space. 


For more options, use:
```bash
raptor update --help
```
For further information and detailed instructions, please consult the comprehensive manual accompanying the HIBF tool.