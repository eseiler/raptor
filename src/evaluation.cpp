#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <regex>
#include <filesystem>

// GLOBAL PARAMETERS
const extern int kmer_size = 20; // For simulated data, use 32, for simulated data, 20
const extern int window_size = kmer_size;
const extern double fpr = 0.05;
const extern int num_hash_functions = 2;
const extern int query_threads = 32;
const extern int query_errors = 2;
const extern bool generate_reads = true;
const extern int number_of_reads = 1000000;
const extern int read_length = 250;
const extern int batch_size = 1;
const extern bool create_queries_once = false;


/*!\brief Reads the time as is stored in the output.time file created by raptor when qureying with the '-time' flag.
 * \param[in] filename
 * \return lastNumber
*/
double read_time(std::string filename){
    std::ifstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cout << "Failed to open the output file!" << std::endl;
    }
    std::string line;
    std::getline(outputFile, line);
    std::getline(outputFile, line);

    size_t lastTabPos = line.find_last_of('\t');
    std::string lastNumberStr = line.substr(lastTabPos + 1);

    double lastNumber{0};     // Convert the extracted substring to a double
    lastNumber = std::stod(lastNumberStr);

    outputFile.close();
    return lastNumber;
}

//!\brief Convert elapsed time string to seconds
double convertElapsedTime(const std::string& elapsedTime) {
    std::regex timeRegex("(\\d+):(\\d+\\.\\d+)");
    std::smatch match;
    if (std::regex_match(elapsedTime, match, timeRegex) && match.size() == 3) {
        int hours = std::stoi(match[0]); // 6:11:09elapsed
        int minutes = std::stoi(match[1]);
        double seconds = std::stod(match[2]);
        return  hours*3600 + minutes * 60 + seconds; // check how time works. hours * 60*60
    } else {
        std::regex timeRegex(R"((\d+):(\d+):(\d+))");
        std::smatch match;

        if (std::regex_match(elapsedTime, match, timeRegex) && match.size() == 4) {
            int hours = std::stoi(match[1]);
            int minutes = std::stoi(match[2]);
            int seconds = std::stoi(match[3]);

            return hours * 3600 + minutes * 60 + seconds;
        }

    }
    return 0.0;
} // TODO make this applicable if runtime takes longer!

//!\brief  Find a word 'SearchName' in the file, and return if it was found.
bool find_rebuild(std::string filename, std::string searchName) {
    bool found = false;
    std::ifstream inputFile(filename);
    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            if (line.find(searchName) != std::string::npos) {
                found = true;
                break;
            }
        }
    }
    return found;
}

//!\brief Executes a command (a string) and stores the output, and extracts time and memory usage.
std::tuple<int, double> execute_command(std::string outputFileName, std::string command) {
    command = "/usr/bin/time " + command + " > " + outputFileName + " 2>&1";
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) { //-f"%M"
        std::cout << "Command execution failed!" << std::endl;
    }
    pclose(pipe);

    std::ifstream outputFile(outputFileName);
    if (!outputFile.is_open()) {
        std::cout << "Failed to open the output file!" << std::endl;
    }
    std::string line;
    std::string new_line;
    std::string previous_line;

    while (std::getline(outputFile, new_line)) {
        line = previous_line;
        previous_line = new_line;
    }
    outputFile.close();

    // Extract elapsed time and maxresident set size using regular expressions
    std::regex timeRegex(".* (\\d+:\\d+\\.\\d+)elapsed.*");//should we change this for long commands that might take over an hour?
    std::regex sizeRegex(".* (\\d+)maxresident.*");
    std::smatch timeMatch, sizeMatch;
    int maxSize; double elapsedSeconds;
    if (std::regex_match(line, timeMatch, timeRegex) && timeMatch.size() == 2) {
        std::string elapsedTime = timeMatch[1];
        elapsedSeconds = convertElapsedTime(elapsedTime);
        std::cout << "Elapsed Time (seconds): " << elapsedSeconds << std::endl;
    }

    if (std::regex_match(line, sizeMatch, sizeRegex) && sizeMatch.size() == 2) {
        maxSize = std::stoi(sizeMatch[1]);
        std::cout << "Maxresident Set Size: " << maxSize << std::endl;
    }else{
        std::cout << "no set size" << std::endl;

    }

    return std::make_tuple(maxSize, elapsedSeconds);
}


//!\brief Write 100 lines from 'filename' to 'filename_queries'.
void write_to_fasta(std::string filename_queries, std::string filename){
    std::ifstream source_file(filename);
    std::ofstream dest_file(filename_queries, std::ios_base::app);
    if (!source_file.is_open() || !dest_file.is_open())  std::cout << "Failed to open the file! " + filename_queries + " or " + filename << std::endl;
    std::string line;
    bool isHeader = false;
    int queryCount=0;
    while (std::getline(source_file, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>')
            continue;
        dest_file << ">" + filename << std::endl;    // write fasta header.
        std::string total_characters;
        dest_file << line << std::endl;
        queryCount +=1;
        if (queryCount == 100){ // the number of queries per file.
            break;
        }
    }
    source_file.close();
    dest_file.close();
    std::cout << "FASTA entries sampled and appended successfully!" << std::endl;
}

//!\brief Writes 'user_bin_filename' to the end of 'existing_filenames'.
void write_to_txt(std::string existing_filenames, std::string user_bin_filename){
    std::ofstream dest_file(existing_filenames, std::ios_base::app);
    if ( !dest_file.is_open())  std::cout << std::endl << "Failed to open the file! " + existing_filenames;
    std::string line;
    dest_file << user_bin_filename << std::endl ;    // write filename TODO it is very important that there are no empty lines, this will cause an error in the layout method.
    dest_file.close();
}

//!\brief Create the queries.
std::tuple<int, std::vector<std::string>> create_query_file(std::string all_paths, std::string filename_queries,
                      std::string folder, std::string tmp_folder = "tmp"){
    int result_delete = std::remove(filename_queries.c_str()); // Delete the existing file, if it is exists
    if (generate_reads){
        std::string filename_executable = folder + "generate_reads";
        std::string command = filename_executable +
                                " --output " + filename_queries +
                                " --errors " +  std::to_string(query_errors) +
                                " --number_of_reads "  +  std::to_string(number_of_reads)  +
                                " --read_length "  +  std::to_string(read_length)  +
                                " " + all_paths;
        std::string ouptut_file = folder + "/evaluation/"+ tmp_folder +"/create_query_file.txt";
        execute_command(ouptut_file, command);
        if (not find_rebuild(ouptut_file, "[SUCCESS]")){
            std::cout << "[ERROR] error message detected!" <<std::flush;
            std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
        }
        std::vector<std::string> dummy_vec = {"dummy string", " "};

        return std::make_tuple(1, dummy_vec);
    }else{
        // 1.1 load all the bin paths
        std::vector<std::string> filenames_existing_bins;      // Vector to store the filenames
        std::ifstream file(all_paths);
        if (!file.is_open()) std::cout << "Failed to open the file! " + all_paths << std::endl;
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty()) {
                filenames_existing_bins.push_back(line);
            }
        }

        // 1.2 create the queries
        std::vector<std::string> tmp_query_filenames;
        for (auto user_bin_filename : filenames_existing_bins){
            write_to_fasta(filename_queries, user_bin_filename);
            // also write single files, such that they can be queried seperately
            size_t lastSlashPos = user_bin_filename.find_last_of("/");
            std::string lastPart = user_bin_filename.substr(lastSlashPos + 1);
            std::string tmp_query_filename = folder + "evaluation/tmp/" + lastPart; // these 3 lines are not important actually.
            write_to_fasta(tmp_query_filename, user_bin_filename);
            tmp_query_filenames.push_back(tmp_query_filename);
        }
            return std::make_tuple((int) filenames_existing_bins.size(), tmp_query_filenames);
    }
}

//!\brief Insert single bin and measure insertion time.
std::tuple<int, double, bool, bool>  insert_ub(std::string filename_ub, std::string filename_index,
                                   std::string filename_executable, std::string sketch_directory,
                                   std::string folder, std::string insertion_method, std::string tmp_folder="tmp"){
    std::string command = filename_executable + // Create the command.
                                    " update " +
                                    " --hibf --insert-UBs " +
                                    " --input " + filename_index +
                                    " --bins " + filename_ub +
                                    " --sketch-directory " + sketch_directory +
                                    " --insertion-method " + insertion_method +
                                    " --sequence-similarity " +
                                    " --output " +  filename_index;
    std::string ouptut_file = folder + "evaluation/" + tmp_folder+"/" + "insertion_output.txt"; // store any output that would otherwise be written to the terminal.
    auto memory_time = execute_command(ouptut_file, command);
    if (find_rebuild(ouptut_file, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::flush;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;
    }
    if (not find_rebuild(ouptut_file, "[SUCCESS]")){
        std::cout << "[ERROR] error message detected!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;
    }
    return std::make_tuple(std::get<0>(memory_time), std::get<1>(memory_time), find_rebuild(ouptut_file, "Full rebuild"), find_rebuild(ouptut_file, "Partial rebuild"));
}

//!\brief Rebuild the index in the naive way.
std::tuple<int, double, bool, bool>  rebuild_index(std::string filename_ub, std::string filename_index,
                                         std::string filename_executable, std::string sketch_directory,
                                         std::string folder, std::string existing_filenames_building, std::string tmp_folder ="tmp"){
    system(("rm " + filename_index).c_str());
    std::string layout_file = folder + "evaluation/" + tmp_folder + "/temporary_layout.txt";
    std::string filename_executable_chopper = folder + "chopper"; //"/mnt/c/Users/myrth/Desktop/coding/raptor/lib/chopper/build/bin/chopper";//folder +  "chopper"; /
    std::string command_build = filename_executable + " build --fpr " +std::to_string(fpr) +
                                                      " --kmer "  + std::to_string(kmer_size) +
                                                      " --window " + std::to_string(window_size) +
                                                      " --hibf --output " + filename_index + " " + layout_file;
    std::string command_layout = filename_executable_chopper + " --num-hash-functions " + std::to_string(num_hash_functions) +
                                                               " --false-positive-rate " + std::to_string(fpr) +
                                                               " --input-file " + existing_filenames_building +
                                                               " --output-filename " + layout_file +
                                                               " --update-UBs 0 " + // no empty bins should be included in the naive method.
                                                               " --output-sketches-to "  + folder + "chopper_sketch_sketches " +
                                                               " --kmer-size " +  std::to_string(kmer_size); //--tmax 64

    std::string ouptut_file = folder + "evaluation/" + tmp_folder + "/rebuild_output.txt";
    auto memory_time_layout = execute_command(ouptut_file, command_layout);
    if (find_rebuild(ouptut_file, "rror")){
        std::cout << "error message!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();

    }
    auto memory_time_build = execute_command(ouptut_file, command_build);
    if (find_rebuild(ouptut_file, "rror")){
        std::cout << "error message!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();

    }
    return std::make_tuple(std::max(std::get<0>(memory_time_layout), std::get<0>(memory_time_build)), // take the maximum when it comes to max residence time.
            std::get<1>(memory_time_layout) + std::get<1>(memory_time_build), 1,0);
}

//!\brief Query multiple samples
std::tuple<int, double>  query_all_ubs(std::string filename_queries, std::string filename_index,
                                       std::string filename_executable,  int number_of_queries,
                                       std::string folder,
                                       std::vector<std::string> tmp_query_filenames,
                                       bool all_at_once = true, std::string tmp_folder="tmp" ){
        std::string filename_ouptut = folder + "evaluation/" + tmp_folder + "/" + "query_result.txt"; // temporary
        std::string command = filename_executable + // Command for querying
                              " search " +
                              " --hibf " +
                              " --error " +  std::to_string(query_errors) +
                              " --index " + filename_index +
                              " --fpr "  +  std::to_string(fpr) +
                              " --query " + filename_queries + // must be a fastq (or a fasta?)  file with all multiple sequences.
                              " --threads " +  std::to_string(query_threads) +
                              " --time --output " + filename_ouptut;// include threads if possible?
        auto memory_time = execute_command(folder + "evaluation/"+ tmp_folder + "/" + "query_output.txt", command);
        if (find_rebuild(filename_ouptut, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::flush;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
        }
        std::string filename_ouptut_time = filename_ouptut + ".time";
        double memory = std::get<0>(memory_time); // This is a maximum, and does not need to be devided by the number of queries. It is also fine that it includes loading the index.
        double time = read_time(filename_ouptut_time); // Extract time from the output file from raptor.
        if (not generate_reads) time = time / number_of_queries; // average. But for generate reads we have always a constant number of queries.
        return std::make_tuple(memory, time);
}


//!\brief Extract filenames from a file with multple paths
std::vector<std::string> extract_filenames(std::string insertion_paths){
    std::ifstream inputFile(insertion_paths);
    if (!inputFile.is_open()) {
        std::cout << "Failed to open the insertion paths file!" << std::endl;
    }
    std::vector<std::string> filenames;
    std::string user_bin_filename;
    while (std::getline(inputFile, user_bin_filename)) {
        filenames.push_back(user_bin_filename);
    }
    inputFile.close();
    return filenames;
}

//!\brief Measure the size of the uncompressed index in kB using a bash code.
int file_size_bash(std::string filename_index){
        system(("du -s " + filename_index + " | cut -f1 > temp_file.txt" ).c_str()); // Stores the file size to a temporary file.
        std::string result;
        std::ifstream temp_file("temp_file.txt");
        if (temp_file) {
            std::getline(temp_file, result);
            temp_file.close();
        }
        if (!result.empty() && result.back() == '\n') result.pop_back(); // Remove newline characters from the end of the result string
        std::remove("temp_file.txt");         // Remove the temporary file
        try{
            int fileSize = std::stod(result); // convert the string to an integer.
            return fileSize;
        }
            catch (const std::exception& e) {
            return 0;
        }
}


//!\brief Create a string of an array to save the results.
template <typename T>
std::string outstring(std::vector<T> out_array){
    std::string out = "";
    for (auto &item: out_array) { std::string x = std::to_string(item); out = out + x + ", "; }
    return out;
}

// SAVE RESULTS
// Saving results to a python library, that can be read in to make figures.
int write_to_python(std::string result_folder,
                    std::vector<double> time_insertion, std::vector<double> time_query,
                    std::vector<double> memory_insertion, std::vector<double> memory_query,
                    std::vector<int> size_index, std::vector<int> size_index_bash, std::vector<int> rebuild){
    std::ofstream o(result_folder + ".py"); // "here.py");//
    o << "time_insertion = [" << outstring(time_insertion) << "]" <<std::endl;
    o << "time_query = [" << outstring(time_query) << "]" <<std::endl;
    o << "memory_insertion = [" << outstring(memory_insertion) << "]" <<std::endl;
    o << "memory_query = [" << outstring(memory_query) << "]" <<std::endl;
    o << "size_index = [" << outstring(size_index) << "]" <<std::endl;
    o << "size_index_bash = [" << outstring(size_index_bash) << "]" <<std::endl;
    o << "rebuild = [" << outstring(rebuild) << "]" <<std::endl;

    o.close();
    return 0;
}



//////////////////////////////////////////////////

int main_insert_ub(){
    // ADDITIONAL INPUT PARAMETERS
    std::string folder = std::filesystem::current_path(); folder += "/";
    std::cout << "Give the test_folder, should be located within the evaluation folder: ";
    std::string input_test; std::cin >> input_test;
    std::cout << "Give the tmp folder, should be located within the evaluation folder: ";
    std::string tmp_folder; std::cin >> tmp_folder;
    std::cout << "Give the insertion method \"find_ibf_size_splitting\", \"find_ibf_idx_traverse_by_fpr\", \"naive\", \"find_ibf_idx_traverse_by_similarity\",  \"find_ibf_idx_ibf_size\",  : ";
    std::string insertion_method; std::cin >> insertion_method;
    std::vector<std::string> insertion_methods = {"find_ibf_size_splitting", "find_ibf_idx_traverse_by_fpr", "naive", "find_ibf_idx_traverse_by_similarity", "find_ibf_idx_ibf_size"};
    if (not (std::find(insertion_methods.begin(), insertion_methods.end(), insertion_method) != insertion_methods.end())) // Check if a correct
        std::cout << "Error: Invalid insertion method.\n";

    // CREATE OTHER VARIABLES.
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt"; // This file contains the paths of the samples in the current index.
    std::string filename_index_original = "hibf.index"; // this could best be an index without empty bins.
    std::string filename_executable = folder +  "raptor"; // Run the code from the same folder where the raptor executable is located.
    std::string result_folder = folder + "evaluation/" + input_test + "/results/" + tmp_folder + "/"; // Here goes the result output, such as the python files.
    std::string sketch_directory = folder + "evaluation/" + tmp_folder  + "/chopper_sketch_sketches";
    std::string filename_queries_existing = folder + "evaluation/" + tmp_folder + "/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/" + tmp_folder + "/" + "queries.fasta";
    std::string existing_filenames_building = folder + "evaluation/" + tmp_folder + "/" + "existing_filenames.txt";

    // MAKE DIRECTORIES, COPY FILES
    system(("mkdir " + tmp_folder + "").c_str()); //
    system(("mkdir -p " + result_folder).c_str());
    system(("mkdir " +folder + "evaluation/" + tmp_folder).c_str());
    std::remove(existing_filenames_building.c_str()); // delete file first
    system(("yes | cp -f " + all_paths + " " + existing_filenames_building).c_str()); //make a copy of the file
    system(("yes | cp -f -r " + folder  + "chopper_sketch_sketches" + " " + sketch_directory).c_str());
    system("cp ../_deps/raptor_chopper_project-src/build/bin/chopper chopper"); // Copy over chopper executable on my own system

    // WRITE VARIABLES TO A CONFIG FILE
    std::ofstream config(folder + "evaluation/" + input_test + "/results/" + tmp_folder + "/" + "config.txt");
    config  << "insertion_method: " << insertion_method <<std::endl; //TEMPORARY
    config  << "folder: " << folder <<std::endl;
    config << "tmp_folder: " << tmp_folder <<std::endl;
    config << "kmer_size: " << kmer_size <<std::endl;
    config << "window_size: " << window_size <<std::endl;
    config << "fpr: " << fpr <<std::endl;
    config << "num_hash_functions: " << num_hash_functions <<std::endl;
    config << "query_threads: " << query_threads <<std::endl;
    config << "query_errors: " << query_errors <<std::endl;
    config << "generate_reads: " << generate_reads <<std::endl;
    config << "number_of_reads: " << number_of_reads <<std::endl;
    config << "read_length: " << read_length <<std::endl<<std::endl;
    config << "insertion_paths: " << insertion_paths <<std::endl;
    config << "all_paths " << all_paths <<std::endl;
    config << "filename_index_original: "<< filename_index_original <<std::endl;
    config << "existing_filenames_building:" << existing_filenames_building <<std::endl;
    config << "filename_queries_existing: " << filename_queries_existing <<std::endl;
    config << "filename_queries: " << filename_queries <<std::endl;
    config << "sketch_directory: " << sketch_directory <<std::endl;
    config << "result_folder: " << result_folder <<std::endl;
    config.close();

    // CREATE QUERIES, EXTRACT FILES.
    auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
    std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);


    // LOOP OVER INSERTION METHODS.     Comment out the loop to only use the specified insertion method.
for (auto insertion_method: {"find_ibf_size_splitting", "find_ibf_idx_traverse_by_fpr", "naive", "find_ibf_idx_traverse_by_similarity",  "find_ibf_idx_ibf_size",   }){
    std::string filename_index =  folder + "evaluation/" + tmp_folder + "/" +insertion_method+ "_" +  filename_index_original;
    int number_of_files = extract_filenames(all_paths).size(); // Count the number of existing files.
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); // make a copy of the file
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); // make a copy of the file
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query; // Create vectors to store the results in
    std::vector<int> size_index, size_index_bash, rebuilds; // Result vectors

    // check time and size before updating
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));
    size_index_bash.push_back(file_size_bash(filename_index));
    int counter=0;
    for (const std::string& user_bin_filename : user_bin_filenames) {
        counter += 1;
        bool query_this_iteration = (counter%batch_size==0);
        number_of_files += 1;
        std::cout << "Method: " << insertion_method << ", with current bin number: " << number_of_files << std::endl << std::flush;
        size_t lastSlashPos = user_bin_filename.find_last_of("/");         // store the filename to insert to a temporary file, 'tmp_filename'.
        std::string lastPart = user_bin_filename.substr(lastSlashPos + 1);
        std::string tmp_filename = folder + "evaluation/" + tmp_folder + "/" + lastPart;
        std::ofstream tmp_insertion_filename(tmp_filename);
        tmp_insertion_filename << user_bin_filename;
        tmp_insertion_filename.close();

        if (not create_queries_once){
            if (generate_reads){         // update all paths and queries
                if (query_this_iteration){
                    auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
             }
            }else{
                write_to_fasta(filename_queries, user_bin_filename); // all queries together
                std::string tmp_query_filename = folder + "evaluation/" + tmp_folder + "/" + lastPart;
                tmp_query_filenames.push_back(tmp_query_filename);
            }
}
        std::tuple<int, double, bool, bool> memory_time_insertion; // tuple to store the results in.
        if (insertion_method != "naive"){
            memory_time_insertion = insert_ub(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, insertion_method, tmp_folder); //  measure insertion time and memory
        }else{
            write_to_txt(existing_filenames_building, user_bin_filename);
            memory_time_insertion = rebuild_index(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, existing_filenames_building, tmp_folder); //  measure insertion time and memory
        }
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        if (std::get<3>(memory_time_insertion) and not std::get<2>(memory_time_insertion)){ // the rebuild array keeps track of partial and full rebuilds occuring.
            rebuilds.push_back(2); // a '2' inidicates a partial rebuild.
        }else{
            rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        }

        if (query_this_iteration){
            auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
            memory_query.push_back(std::get<0>(memory_time_queries));
            time_query.push_back(std::get<1>(memory_time_queries));
            size_index_bash.push_back(file_size_bash(filename_index));
        }
        if (counter%50==0){ // Store intermediate results once in the 50 insertions.
            std::cout << "Saving intermediate results. " <<std::endl;
            write_to_python(result_folder + insertion_method, time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
        }

    }
    if (batch_size > 1){ // make sure you also query the index at the end, if you use batches.
        if (not create_queries_once){
            if (generate_reads){
                auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
            }
                auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
                memory_query.push_back(std::get<0>(memory_time_queries));
                time_query.push_back(std::get<1>(memory_time_queries));
                size_index_bash.push_back(file_size_bash(filename_index));
        }
    }

    write_to_python(result_folder + insertion_method, time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
}
    return 0;
}

//INSERT SEQUENCES
//1. Insert single bin and measure insertion time.
std::tuple<int, double, bool, bool>  insert_sequences(std::string filename_ub, std::string filename_index,
                                   std::string filename_executable, std::string sketch_directory,
                                   std::string folder, std::string tmp_folder = "tmp"){
    std::string command = filename_executable +
                            " update " +
                        " --hibf --insert-sequences" +
                        " --input " +
                          filename_index +
                          " --bins " +
                            filename_ub +
                            " --sketch-directory " +
                            sketch_directory +
                            " --sequence-similarity " +
                            " --empty-bin-sampling 0.1 " + // some empty bins to aid partial rebuilds.
                          " --output " +
                          filename_index;
    //--insert_sequence_appendix. By default the ending is "_insertsequences". It is the responsibility of the user to update the fasta files themselves, such as by using the concatenate cat command.
    //  a list of their filenames through the paremeter --bins

    //1. Build a large index, with both (?) 64 and 1024 bins -> store this index in a folder "mixed_bins".
    // (later, bins of different sizes). Then add all sequence content of bin_00 to bin_0000, line by line (or some other bigger file). Concatenate the lines of bin_00 to bin_0000
    //2. Add to a larger bin e.g. bin_00.
    //"bin_0000_insertsequences"

    std::string insert_sequence_appendix = "_insertsequences";
    std::string ouptut_file = folder + "evaluation/" + tmp_folder + "/" + "insertion_output.txt";
    auto memory_time = execute_command(ouptut_file, command);
    if (find_rebuild(ouptut_file, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::endl;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    if (not find_rebuild(ouptut_file, "[SUCCESS]")){
        std::cout << "[ERROR] error message detected!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    return std::make_tuple(std::get<0>(memory_time), std::get<1>(memory_time), find_rebuild(ouptut_file, "Full rebuild"), find_rebuild(ouptut_file, "Partial rebuild"));
}
//2. Insert single bin and measure insertion time.


int main_insert_seq(){
    std::cout << "WARNING: Make sure you set the bin_0000 back to its original state ";
    std::string folder = std::filesystem::current_path(); folder += "/";
    std::cout << "Give the test_folder, should be located within the evaluation folder: . Also 'bin_bonus.fasta' should be in the test folder. ";
    std::string input_test; std::cin >> input_test;
    std::cout << "Give the tmp folder: ";
    std::string tmp_folder; std::cin >> tmp_folder;
    int number_of_lines  = 95872; // number of lines at a time.
// PARAMETERS
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";"evaluation.index"; // this could best be an index without empty bins.

    std::string filename_executable = folder +  "raptor";

    // output
    system(("mkdir " + folder + "evaluation/" + input_test + "/results/").c_str());
    std::string result_folder = folder + "evaluation/" + input_test + "/results/results";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";
    system("mkdir evaluation");
    system("mkdir evaluation/results");
    std::string filename_queries_existing = folder + "evaluation/"+ tmp_folder +"/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/tmp/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/"+ tmp_folder ).c_str());
    std::string existing_filenames_building = folder + "evaluation/"+ tmp_folder +"/" + "existing_filenames.txt";
    std::remove(existing_filenames_building.c_str()); // delete file first
    system(("yes | cp -f " + all_paths + " " + existing_filenames_building).c_str()); //make a copy of the file



    /////////////
    //if (isLastLineEmpty(all_paths))   deleteLastLine(all_paths);
    //copy chopper
    system("cp ../_deps/raptor_chopper_project-src/build/bin/chopper chopper");

    std::cout << folder <<std::flush;
    std::cout << insertion_paths <<std::flush;
    std::cout << all_paths <<std::flush;
    std::cout << filename_index_original <<std::flush;
//    std::cout << filename_index <<std::flush;
    std::cout << existing_filenames_building <<std::flush;
    std::cout << filename_queries_existing <<std::flush;
    std::cout << filename_queries <<std::flush;
    std::cout << sketch_directory <<std::flush;
    std::cout << result_folder <<std::flush;



    auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
    int number_of_files = extract_filenames(existing_filenames_building).size(); std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::cout << "number of lines: " << number_of_files << std::endl;
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);

    std::string filename_index =  folder + "evaluation/"+ tmp_folder +"/" +"sequence_insertions"+ "_" +  filename_index_original;

    std::cout << std::endl << std::endl  << "sequence_insertions" << std::endl  << std::flush;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file

    // Result vectors
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds;
    // check time and size before updating
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));
    size_index_bash.push_back(file_size_bash(filename_index));

        system("mkdir tmp"); //


        //create file with sequence content to insert "bin_0000_insertsequences"
        // pick next 10 lines from bin_0000

    std::string file_to_insert_from = folder + "evaluation/" + input_test + "/bin_bonus.fasta";
//    std::ifstream inputFile(file_to_insert_from);  // Open input file
//    if (!inputFile.is_open()) {
//        std::cerr << "Failed to open input file." << std::endl;
//        return 1;
//    }
    std::ifstream opened_file(all_paths);
    std::string insert_to_ub;
    std::getline(opened_file, insert_to_ub);
    std::cout << "the sample to insert to: " << insert_to_ub;
    // WHEN TESTING ON REAL DATA: 1. copy all_paths to the tmp folder (all_paths_original and all_paths). 2. change the first line. 3. copy over the sample in the first line.
    // sample 10% of the chosen ub.

    std::string insert_to_ub_paths = "evaluation/"+ tmp_folder + "/insert_to_ub_paths.txt";
    std::ofstream pathFile(insert_to_ub_paths);
    pathFile << insert_to_ub << std::endl;
    pathFile.close();
    std::string filename_new_sequences = insert_to_ub.substr(0, insert_to_ub.find_last_of('.')) + "_insertsequences"
            + insert_to_ub.substr(insert_to_ub.find_last_of('.'));

    int fileCount = 0;
    int lineCount = 1;
    std::string line;

    std::ofstream outputFile(filename_new_sequences, std::ios_base::trunc);  // Open output file in append mode
    outputFile << ">Header1\n";
    outputFile.close();
    std::system(("cp " + insert_to_ub + " " +  insert_to_ub + "_original").c_str() ); // if only 1 file is updated.

    for (int _ = 1; _ <= 500; ++_){
        std::system(("sed -n '" + std::to_string(lineCount) + ","+ std::to_string(lineCount + number_of_lines)
        + "p' " + file_to_insert_from + " >> " + filename_new_sequences).c_str() ); // if only 1 file is updated.
//        std::system(("[[ -n \""+filename_new_sequences+"\" && \"$(tail -n 1 \""+filename_new_sequences+"\")\" == '>'* ]] && sed -i '$d' \""+filename_new_sequences+"\"").c_str()); // removes last line if it starts with ">"
//        std::system(("[[ -n \""+filename_new_sequences+"\" && \"$(sed '2q;d' \""+filename_new_sequences+"\")\" == '>'* ]] && sed -i '2d' \""+filename_new_sequences+"\"").c_str()); // removes second line if it starts with ">"
        //std::system(("echo -e \"\\n\" >> " + filename_new_sequences).c_str()); // removes second line if it starts with ">"
        std::system(("echo \"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\" >> " + filename_new_sequences).c_str()); // adds a small dummy sequence to prevent the file from ending with a header-line.


        lineCount += number_of_lines;

//    while (getline(inputFile, line)) { // TODO a bash script would be faster to copy over the lines x to x + number_of_lines. Where x is 0 and increases gradually.
//        if (line.empty() || line[0] == '>') {
//            continue;  // Skip empty lines or lines starting with '>'
//        }
//        if (lineCount % number_of_lines == 0 and lineCount!= 0) { // 10 lines will be inserted a time.
        std::cout << "lineCount: " << lineCount << ", insertion #" << lineCount/number_of_lines << std::endl;
        std::system(("cat " + filename_new_sequences + " >> " + insert_to_ub ).c_str());


         std::tuple<int, double, bool, bool> memory_time_insertion = insert_sequences(insert_to_ub_paths, filename_index,
                               filename_executable, sketch_directory,
                               folder, tmp_folder);
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        if (std::get<3>(memory_time_insertion) and not std::get<2>(memory_time_insertion)){
            rebuilds.push_back(2); // 2 inidicates a partial rebuild.
        }else{
            rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        }

        auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));
        size_index_bash.push_back(file_size_bash(filename_index));

        std::remove(filename_new_sequences.c_str());         // empty the file outputFileName, such that 10 new sequences can be inserted.
        std::ofstream outputFile(filename_new_sequences);  // Open output file in append mode
        outputFile << ">Header\n";     outputFile.close();



        if (lineCount%(number_of_lines*20)==0){
            std::cout << "saving intermediate results" <<std::endl;
                    write_to_python(result_folder + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
        }

//        std::ofstream outputFile(filename_new_sequences, std::ios_base::app);  // Open output file in append mode
//        outputFile << line << std::endl; outputFile.close();
//        lineCount++;
    }

    write_to_python(result_folder + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
    std::system(("mv " + insert_to_ub + "_original " + insert_to_ub + " ").c_str() );

    return 0;
    //system
    // at the end set back the original file.
    // tmp copy

}



int main_insert_seq_2(){ // insert sequence material from seperate samples, one by one .
    std::cout << "WARNING: Make sure you set the bin_0000 back to its original state ";
    std::string folder = std::filesystem::current_path(); folder += "/";
    std::cout << "Give the test_folder, should be located within the evaluation folder: . Also 'bin_bonus.fasta' should be in the test folder. ";
    std::string input_test; std::cin >> input_test;
    std::cout << "Give the tmp folder: ";
    std::string tmp_folder; std::cin >> tmp_folder;
    int number_of_lines  = 95872; // number of lines at a time.
// PARAMETERS
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";"evaluation.index"; // this could best be an index without empty bins.

    std::string filename_executable = folder +  "raptor";

    // output
    system(("mkdir " + folder + "evaluation/" + input_test + "/results/").c_str());
    std::string result_folder = folder + "evaluation/" + input_test + "/results/results";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";

    std::string filename_queries_existing = folder + "evaluation/"+ tmp_folder +"/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/tmp/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/"+ tmp_folder ).c_str());
    std::string existing_filenames_building = folder + "evaluation/"+ tmp_folder +"/" + "existing_filenames.txt";
    std::remove(existing_filenames_building.c_str()); // delete file first
    system(("yes | cp -f " + all_paths + " " + existing_filenames_building).c_str()); //make a copy of the file



    /////////////
    //if (isLastLineEmpty(all_paths))   deleteLastLine(all_paths);
    //copy chopper
    system("cp ../_deps/raptor_chopper_project-src/build/bin/chopper chopper");

    std::cout << folder <<std::flush;
    std::cout << insertion_paths <<std::flush;
    std::cout << all_paths <<std::flush;
    std::cout << filename_index_original <<std::flush;
//    std::cout << filename_index <<std::flush;
    std::cout << existing_filenames_building <<std::flush;
    std::cout << filename_queries_existing <<std::flush;
    std::cout << filename_queries <<std::flush;
    std::cout << sketch_directory <<std::flush;
    std::cout << result_folder <<std::flush;



    auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
    int number_of_files = extract_filenames(existing_filenames_building).size(); std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::cout << "number of lines: " << number_of_files << std::endl;
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);

    std::string filename_index =  folder + "evaluation/"+ tmp_folder +"/" +"sequence_insertions"+ "_" +  filename_index_original;

    std::cout << std::endl << std::endl  << "sequence_insertions" << std::endl  << std::flush;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file

    // Result vectors
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds;
    // check time and size before updating
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));
    size_index_bash.push_back(file_size_bash(filename_index));

    std::ifstream opened_file(all_paths);
    std::string insert_to_ub;
    std::getline(opened_file, insert_to_ub);
    std::cout << "the sample to insert to: " << insert_to_ub;
    // WHEN TESTING ON REAL DATA: 1. copy all_paths to the tmp folder (all_paths_original and all_paths). 2. change the first line. 3. copy over the sample in the first line.
    // sample 10% of the chosen ub.

    std::string insert_to_ub_paths = "evaluation/"+ tmp_folder + "/insert_to_ub_paths.txt";
    std::ofstream pathFile(insert_to_ub_paths);
    pathFile << insert_to_ub << std::endl;
    pathFile.close();
    std::string filename_new_sequences = insert_to_ub.substr(0, insert_to_ub.find_last_of('.')) + "_insertsequences"
            + insert_to_ub.substr(insert_to_ub.find_last_of('.'));

    int fileCount = 0;
//
//    std::ofstream outputFile(filename_new_sequences, std::ios_base::trunc);  // Open output file in append mode
//    outputFile << ">Header1\n";
//    outputFile.close();
//    std::system(("cp " + insert_to_ub + " " +  insert_to_ub + "_original").c_str() ); // if only 1 file is updated.


    for (const std::string& user_bin_filename : user_bin_filenames) {
        fileCount += 1;

//        if (line.empty() || line[0] == '>') {
//            continue;  // Skip empty lines or lines starting with '>'
//        }
        //if (lineCount % number_of_lines == 0 and lineCount!= 0) { // 10 lines will be inserted a time.
        std::cout << "fileCount: " << fileCount << std::endl;
        std::system(("cp " + user_bin_filename + " " + filename_new_sequences ).c_str());

        std::system(("cat " + filename_new_sequences + " >> " + insert_to_ub ).c_str());


         std::tuple<int, double, bool, bool> memory_time_insertion = insert_sequences(insert_to_ub_paths, filename_index,
                               filename_executable, sketch_directory,
                               folder, tmp_folder);
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        if (std::get<3>(memory_time_insertion) and not std::get<2>(memory_time_insertion)){
            rebuilds.push_back(2); // 2 inidicates a partial rebuild.
        }else{
            rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        }

        auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));
        size_index_bash.push_back(file_size_bash(filename_index));



        if (fileCount%(20)==0){
            std::cout << "saving intermediate results" <<std::endl;
                    write_to_python(result_folder + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
        }

    }


    write_to_python(result_folder + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
    std::system(("mv " + insert_to_ub + "_original " + insert_to_ub + " ").c_str() );

    return 0;
    //system
    // at the end set back the original file.
    // tmp copy

}


//DELETE UBs
//1. Insert single bin and measure insertion time.
std::tuple<int, double, bool> delete_ub(std::string filename_ub, std::string filename_index,
                                   std::string filename_executable, std::string sketch_directory,
                                   std::string folder,  std::string tmp_folder){
    std::string command = filename_executable +
                            " update " +
                        " --hibf --delete-UBs" +
                        " --input " +
                          filename_index +
                          " --bins " +
                            filename_ub +
                            " --sequence-similarity " +
                            " --empty-bin-sampling 0.1 " + // some empty bins to aid partial rebuilds.
                          " --output " +
                          filename_index;

    std::string insert_sequence_appendix = "_insertsequences";
    std::string ouptut_file = folder + "evaluation/" + tmp_folder + "/" + "deletion_output.txt";
    auto memory_time = execute_command(ouptut_file, command);
    if (find_rebuild(ouptut_file, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::flush;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    if (not find_rebuild(ouptut_file, "[SUCCESS]")){
        std::cout << "[ERROR] error message detected!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    return std::make_tuple(std::get<0>(memory_time), std::get<1>(memory_time), find_rebuild(ouptut_file, "Svenja+Myrthe"));
}

int main_del_seq(){
    std::cout<< "del ";
    std::string folder = std::filesystem::current_path(); folder += "/";
    std::cout << "Give the test_folder, should be located within the evaluation folder: ";
    std::string input_test; std::cin >> input_test;
    std::cout << "Give the tmp folder: ";
    std::string tmp_folder; std::cin >> tmp_folder;
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";// this could best be an index without empty bins. // start with an index or around 100 bins.
    std::string filename_executable = folder +  "raptor";
    std::string result_folder = folder + "evaluation/" + input_test + "/results/";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";
    system("mkdir evaluation");
    system("mkdir evaluation/results");
    system(("mkdir " + tmp_folder + "").c_str()); //
    system(("mkdir -p " + result_folder).c_str());
    std::string filename_queries_existing = folder + "evaluation/" + tmp_folder + "/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/" + tmp_folder + "/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/" + tmp_folder + "").c_str());
    std::string existing_filenames_building = folder + "evaluation/" + tmp_folder + "/" + "existing_filenames.txt";
    std::remove(existing_filenames_building.c_str()); // delete file first
    system(("yes | cp -f " + all_paths + " " + existing_filenames_building).c_str()); //make a copy of the file
    system("cp ../_deps/raptor_chopper_project-src/build/bin/chopper chopper");

    std::ofstream config(folder + "evaluation/" + input_test + "/results/" + tmp_folder + "/" + "config.txt");
    config  << "folder: " << folder <<std::endl;
    config << "tmp_folder: " << tmp_folder <<std::endl;
    config << "kmer_size: " << kmer_size <<std::endl;
    config << "window_size: " << window_size <<std::endl;
    config << "fpr: " << fpr <<std::endl;
    config << "num_hash_functions: " << num_hash_functions <<std::endl;
    config << "query_threads: " << query_threads <<std::endl;
    config << "query_errors: " << query_errors <<std::endl;
    config << "generate_reads: " << generate_reads <<std::endl;
    config << "number_of_reads: " << number_of_reads <<std::endl;
    config << "read_length: " << read_length <<std::endl<<std::endl;
    config << "insertion_paths: " << insertion_paths <<std::endl;
    config << "all_paths " << all_paths <<std::endl;
    config << "filename_index_original: "<< filename_index_original <<std::endl;
    config << "existing_filenames_building:" << existing_filenames_building <<std::endl;
    config << "filename_queries_existing: " << filename_queries_existing <<std::endl;
    config << "filename_queries: " << filename_queries <<std::endl;
    config << "sketch_directory: " << sketch_directory <<std::endl;
    config << "result_folder: " << result_folder <<std::endl;
    config.close();

    auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file

    int number_of_files = extract_filenames(existing_filenames_building).size(); std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);
    std::string tmp_filename = folder + "evaluation/" + tmp_folder + "/tmp_input_file.txt";
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds, del_or_insert;
    std::string insertion_paths_copy = "evaluation/" + tmp_folder + "/insertion_paths_copy.txt ";
    system(("yes | cp -f " + insertion_paths + " " + insertion_paths_copy).c_str()); //make a copy of the file
    std::string filename_index =  folder + "evaluation/" + tmp_folder + "/" "deletions" + "_" +  filename_index_original;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file


    int user_bin_filenames_counter =0;
    for (int number_of_operatons = 0; number_of_operatons < 1000; number_of_operatons++){
        if (user_bin_filenames_counter + number_of_operatons < user_bin_filenames.size()){
            std::cout <<"number of operation in the series: " << number_of_operatons << std::endl;
    for (int _ = 0; _ < number_of_operatons; _++){
        std::cout <<"number of operation in the series: " << number_of_operatons << ", deletion: "<< _ << std::endl;

        std::ifstream existing_paths_opened(existing_filenames_building);
        std::string user_bin_filename;
        std::getline(existing_paths_opened, user_bin_filename);             // Read the first line from the input file
        existing_paths_opened.close();
        std::ofstream tmp_insertion_filename(tmp_filename);
        tmp_insertion_filename << user_bin_filename;
        tmp_insertion_filename.close();

        auto memory_time_insertion = delete_ub(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, tmp_folder); //  measure insertion time and memory
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        del_or_insert.push_back(1);
        system(("tail -n +2 " + existing_filenames_building +" > tmp.txt").c_str()); // Read the content of the file, excluding the first line /std::to_string(number_of_operatons + 1) +
        system(("mv tmp.txt " + existing_filenames_building).c_str()); // Overwrite the original file with the modified content ==>         delete last from existing_paths.
        memory_query.push_back(-100000); // write a dummy result to the array.
        time_query.push_back(-100000);
        size_index_bash.push_back(file_size_bash(filename_index));
    }

    if (generate_reads){
        auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
        }
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
    if (number_of_operatons and memory_query.size()){
        memory_query.pop_back(); time_query.pop_back();
        }
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));

    for (int _ = 0; _ < number_of_operatons; _++){
        std::string user_bin_filename = user_bin_filenames[user_bin_filenames_counter];
         user_bin_filenames_counter +=1;

        // write new bin to temporary file.
         std::ofstream tmp_insertion_filename(tmp_filename);
         tmp_insertion_filename << user_bin_filename;
         tmp_insertion_filename.close();

         // write to existing filenames.
         std::ofstream existing_paths_opened(existing_filenames_building, std::ios_base::app);
         existing_paths_opened << std::endl <<  user_bin_filename ;
         existing_paths_opened.close();

        auto memory_time_insertion = insert_ub(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, "find_ibf_size_splitting", tmp_folder); //  measure insertion time and memory
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        del_or_insert.push_back(0);
        memory_query.push_back(-100000);
        time_query.push_back(-100000);
        size_index_bash.push_back(file_size_bash(filename_index));

}
     if (generate_reads){
            auto result = create_query_file(existing_filenames_building, filename_queries_existing, folder, tmp_folder);
        }
        memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
        if (number_of_operatons and memory_query.size()){
                    memory_query.pop_back(); time_query.pop_back();
        }
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));
     }
     // expected result: rebuild after some time because merged bins grow in FPR.
     // before that, querying might take longer because of false hits. (depends on how much you query and if you insert the UBs that you deleted, how similart they are)
     write_to_python(result_folder + "deletions", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
//TODO also write del_or_inserts to python file.
 }
}


int main(){
    std::cout<< "seq/ub/del/seq_real_data" ;
    std::string test_type; std::cin >> test_type;
    if (test_type == "ub"){
        main_insert_ub();
    }else if (test_type == "seq"){
        main_insert_seq();
    }else if (test_type == "del"){
        main_del_seq();
    }else if (test_type == "seq_real_data"){
        main_insert_seq_2();}
    return 0;
}
