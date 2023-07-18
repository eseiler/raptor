// --------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
// --------------------------------------------------------------------------------------------------


#include <chrono>
#include <time.h>
#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <regex>
#include <cassert>
#include <filesystem>

#include <raptor/update/load_hibf.hpp>

const extern int kmer_size = 20; // for real data
const extern int window_size = kmer_size;
const extern int fpr = 0.05; // this is again hard-coded, because std::to_string(fpr) converts the fpr into 0.
const extern int num_hash_functions = 2;
const extern int query_threads = 32;
const extern int query_errors = 2;
const extern bool generate_reads = true;
const extern int number_of_reads = 1000000;
const extern int read_length = 250;

#include <iostream>
#include <fstream>
#include <string>

//int read_config() {
//    // Read the configuration from the text file
//    std::ifstream ifs("config.txt");
//    std::string line;
//    while (std::getline(ifs, line)) {
//        std::string key, value;
//        std::size_t pos = line.find('=');
//        if (pos != std::string::npos) {
//            key = line.substr(0, pos);
//            value = line.substr(pos + 1);
//        }
//
//        // Process the configuration variables
//        if (key == "kmer_size") {
//            int kmer_size = std::stoi(value);
//            std::cout << "kmer_size: " << kmer_size << std::endl;
//        } else if (key == "window_size") {
//            int window_size = std::stoi(value);
//            std::cout << "window_size: " << window_size << std::endl;
//        } else if (key == "fpr") {
//            double fpr = std::stod(value);
//            std::cout << "fpr: " << fpr << std::endl;
//        } else if (key == "num_hash_functions") {
//            int num_hash_functions = std::stoi(value);
//            std::cout << "num_hash_functions: " << num_hash_functions << std::endl;
//        } else if (key == "query_threads") {
//            int query_threads = std::stoi(value);
//            std::cout << "query_threads: " << query_threads << std::endl;
//        } else if (key == "query_errors") {
//            int query_errors = std::stoi(value);
//            std::cout << "query_errors: " << query_errors << std::endl;
//        } else if (key == "generate_reads") {
//            bool generate_reads = (value == "true");
//            std::cout << "generate_reads: " << std::boolalpha << generate_reads << std::endl;
//        } else if (key == "number_of_reads") {
//            int number_of_reads = std::stoi(value);
//            std::cout << "number_of_reads: " << number_of_reads << std::endl;
//        } else if (key == "read_length") {
//            int read_length = std::stoi(value);
//            std::cout << "read_length: " << read_length << std::endl;
//        }
//    }
//
//    ifs.close();
//    return 0;
//}
//


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

    // Convert the extracted substring to a double
    double lastNumber{0};
    lastNumber = std::stod(lastNumberStr);

    outputFile.close();
    return lastNumber;
}

// Function to convert elapsed time string to seconds
double convertElapsedTime(const std::string& elapsedTime) {
    std::regex timeRegex("(\\d+):(\\d+\\.\\d+)");
    std::smatch match;

    if (std::regex_match(elapsedTime, match, timeRegex) && match.size() == 3) {
        // int hours = std::stoi(match[0]);
        int minutes = std::stoi(match[1]);
        double seconds = std::stod(match[2]);

        return  minutes * 60 + seconds; // check how time works. hours * 60*60
    }

    return 0.0;
}

// Function to convert maxresident set size string to numeric value
int convertMaxResidentSize(const std::string& maxResidentSize) {
    std::regex sizeRegex("(\\d+)maxresident");
    std::smatch match;

    if (std::regex_match(maxResidentSize, match, sizeRegex) && match.size() == 2) {
        return std::stoi(match[1]);
    }

    return 0;
}

bool find_rebuild(std::string filename, std::string searchName) {
       // is print if it rebuilds, and if so with how many user bins/ at what level

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



void write_to_fasta(std::string filename_queries, std::string filename, double sample_percentage = 0.05){

    std::ifstream source_file(filename);
    std::ofstream dest_file(filename_queries, std::ios_base::app);

    if (!source_file.is_open() || !dest_file.is_open())  std::cout << "Failed to open the file! " + filename_queries + " or " + filename << std::endl;
    std::string line;
    bool isHeader = false;
//    int totalSequenceLines = 0;
//    // Count total number of sequence lines (without headers)
//    while (std::getline(source_file, line)) {
//        if (line.empty())
//            continue;
//        if (line[0] == '>') {
//            isHeader = true;
//            continue;
//        }
//        if (isHeader) {
//            isHeader = false;
//            continue;
//        }
//        totalSequenceLines++;
//    }
//    // Calculate target sample size
//    //targetSampleSize = static_cast<int>(totalSequenceLines * sample_percentage);
//
//    // Reset file pointer to the beginning
//    source_file.clear();
//    source_file.seekg(0, std::ios::beg);


    int queryCount=0;
    while (std::getline(source_file, line)) {
        // for 1 in 100
        //for (int split = 0; split < number_of_splits; split++){
        if (line.empty())
            continue;
        if (line[0] == '>')
            continue;
        dest_file << ">" + filename << std::endl;    // write fasta header.
        std::string total_characters;
//        total_characters += line.substr(0, 70);
//            while (std::getline(source_file, line) && total_characters.length() < 250) {
//                total_characters += line.substr(0, 70);
//            }
//            if (total_characters.length() > 250)
//                total_characters = total_characters.substr(0,250);
//        dest_file << total_characters << std::endl; // write a line of 250 bases by writing 3 lines.
        dest_file << line << std::endl;
        queryCount +=1;
        if (queryCount == 100){ // the number of queries per file.
            break;
        }
    }
    // Random number generator for sampling lines
//    std::random_device rd;
//    std::mt19937 generator(rd());
//    std::uniform_real_distribution<double> distribution(0.0, 1.0);
//    while (std::getline(source_file, line)) {
//        if (line.empty())
//            continue;
//        if (line[0] == '>')
//            continue;
//        if (distribution(generator) <= sample_percentage) {       // Sample x% of the lines
//            dest_file << ">" + filename << std::endl;    // write fasta header.
//            dest_file << line << std::endl;
//        }
//        }
    source_file.close();
    dest_file.close();
    std::cout << "FASTA entries sampled and appended successfully!" << std::endl;
}



bool isLastLineEmpty(const std::string& filename) {
    std::ofstream fileOut(filename);
    fileOut << std::endl;
    fileOut.close();

    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string lastLine;
    std::string line;

    while (std::getline(file, line)) {
        if (!line.empty() && !std::all_of(line.begin(), line.end(), [](unsigned char c) { return std::isspace(c); })) {
            lastLine = line;
        }
    }

    file.close();

    return lastLine.empty();
}

void deleteLastLine(const std::string& filename) {
    std::ifstream fileIn(filename);
    if (!fileIn.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    std::string fileContent;
    std::string line;

    while (std::getline(fileIn, line)) {
        fileContent += line + '\n';
    }

    fileIn.close();

    std::ofstream fileOut(filename);
    if (!fileOut.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    if (!fileContent.empty()) {
        fileContent.pop_back(); // Remove the last newline character
        fileOut << fileContent;
    }

    fileOut.close();
}

void write_to_txt(std::string existing_filenames, std::string user_bin_filename){
    std::ofstream dest_file(existing_filenames, std::ios_base::app);

    if ( !dest_file.is_open())  std::cout << std::endl << "Failed to open the file! " + existing_filenames;


    std::string line;
    dest_file << user_bin_filename << std::endl ;    // write filename TODO it is very important that there are no empty lines, this will cause an error in the layout method.
    dest_file.close();
}



std::tuple<int, std::vector<std::string>> create_query_file(std::string all_paths, std::string filename_queries,
                      double sample_percentage, std::string folder, std::string tmp_folder = "tmp"){

    int result_delete = std::remove(filename_queries.c_str()); // Delete the file
    if (result_delete == 0)  printf("File deleted successfully.\n");
    else printf("Failed to delete the file.\n");

    if (generate_reads){
        std::string filename_executable = folder + "generate_reads";
        std::string command = filename_executable +
                                " --output " + filename_queries +
                                " --errors " +  std::to_string(query_errors) +
                                " --number_of_reads "  +  std::to_string(number_of_reads)  +
                                " --read_length "  +  std::to_string(read_length)  +
                                " " + all_paths;
        std::string ouptut_file = folder + "/evaluation/"+ tmp_folder +"/query_output.txt";
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
            write_to_fasta(filename_queries, user_bin_filename, sample_percentage);
            // also write single files, such that they can be queried seperately
            size_t lastSlashPos = user_bin_filename.find_last_of("/");
            std::string lastPart = user_bin_filename.substr(lastSlashPos + 1);

            std::string tmp_query_filename = folder + "evaluation/tmp/" + lastPart; // these 3 lines are not important actually.
            write_to_fasta(tmp_query_filename, user_bin_filename, sample_percentage);
            tmp_query_filenames.push_back(tmp_query_filename);
        }
            return std::make_tuple((int) filenames_existing_bins.size(), tmp_query_filenames);
    }
}

//1. Insert single bin and measure insertion time.
std::tuple<int, double, bool, bool>  insert_ub(std::string filename_ub, std::string filename_index,
                                   std::string filename_executable, std::string sketch_directory,
                                   std::string folder, std::string insertion_method, std::string tmp_folder="tmp"){
    std::string command = filename_executable +
                                    " update " +
                                    " --hibf --insert-UBs " +
                                    " --input " +
                                    filename_index +
                                    " --bins " +
                                    filename_ub +
                                    " --sketch-directory " +
                                    sketch_directory +
                                    " --insertion-method " +
                                    insertion_method +
                                    " --sequence-similarity " +
                                    " --output " +
                                    filename_index;
    std::string ouptut_file = folder + "evaluation/"+ tmp_folder+"/" + "insertion_output.txt";
    auto memory_time = execute_command(ouptut_file, command);
    if (find_rebuild(ouptut_file, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::flush;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    if (not find_rebuild(ouptut_file, "[SUCCESS]")){
        std::cout << "[ERROR] error message detected!" <<std::flush;
        std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
    }
    return std::make_tuple(std::get<0>(memory_time), std::get<1>(memory_time), find_rebuild(ouptut_file, "Svenja+Myrthe"), find_rebuild(ouptut_file, "Partial rebuild"));
}

//1. Insert single bin and measure insertion time.
std::tuple<int, double, bool, bool>  rebuild_index(std::string filename_ub, std::string filename_index,
                                         std::string filename_executable, std::string sketch_directory,
                                         std::string folder, std::string existing_filenames_building, std::string tmp_folder ="tmp"){
    system(("rm " + filename_index).c_str());
    std::string layout_file = folder + "evaluation/" + tmp_folder + "/temporary_layout.txt";
    std::string filename_executable_chopper = folder + "chopper"; //"/mnt/c/Users/myrth/Desktop/coding/raptor/lib/chopper/build/bin/chopper";//folder +  "chopper"; /
    std::string command_build = filename_executable + " build --fpr " + "0.05" +
                                                      " --kmer "  + std::to_string(kmer_size) +
                                                      " --window " + std::to_string(window_size) +
                                                      " --hibf --output " + filename_index + " " + layout_file;
    std::string command_layout = filename_executable_chopper + " --num-hash-functions " + std::to_string(num_hash_functions) +
                                                               " --false-positive-rate " + "0.05" +
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

//!\brief query a single bin.
std::tuple<int, double>  query_all_ubs(std::string filename_queries, std::string filename_index,
                                       std::string filename_executable,  int number_of_queries,
                                       std::string folder,
                                       std::vector<std::string> tmp_query_filenames,
                                       bool all_at_once = true, std::string tmp_folder="tmp" ){
    if (all_at_once) {
        std::string filename_ouptut = folder + "evaluation/" + tmp_folder + "/" + "query_result.txt"; // temporary
        std::string command = filename_executable +
                              " search " +
                              " --hibf " +
                              " --error " +  std::to_string(query_errors) +
                              " --index " +
                              filename_index +
                              " --fpr "  +  "0.05" +
                              " --query " +
                              filename_queries +
                              // must be a fastq (or a fasta?)  file with all multiple sequences.v
                              " --threads " +  std::to_string(query_threads) +
                              " --time --output "  + // include threads if possible?
                              filename_ouptut;
        auto memory_time = execute_command(folder + "evaluation/"+ tmp_folder + "/" + "query_output.txt", command);
        if (find_rebuild(filename_ouptut, "rror")){
                std::cout << "[ERROR] error message detected!" <<std::flush;
                std::cin.clear(); std::cin.get(); int n; std::cin >> n;//std::exit();
        }
        std::string filename_ouptut_time = filename_ouptut + ".time";
        double memory = std::get<0>(memory_time); // This is a maximum, and does not need to be devided by the number of queries. It is also fine that it includes loading the index.
        //double time = std::get<1>(memory_time)/number_of_queries; // Extract time and memory consupmtion
        double time = read_time(filename_ouptut_time); // Instead Extract time from the output file from raptor.
        if (not generate_reads) time = time / number_of_queries; // average. But for generate reads we have always a constant number of queries.
        return std::make_tuple(memory, time);
    }
    else{
        double mean = 0;
        auto variance_func = [&mean, &number_of_queries](double accumulator, const double & val) {
            return accumulator + ((val - mean)*(val - mean) / (number_of_queries - 1));
        };

        std::vector<double> memories, times;
        for (auto & tmp_query_filename : tmp_query_filenames){
            std::string filename_ouptut = folder + "evaluation/tmp/" + "query_result.txt"; // temporary
            std::string command = filename_executable +
                                  " search " +
                                  " --hibf " +
                                  " --error 1 " +
                                  " --index " +
                                  filename_index +
                                  " --fpr 0.05" +
                                  " --query " +
                                    tmp_query_filename + // must be a fastq (or a fasta?)  file with all multiple sequences.v
                                  " --output " +
                                  filename_ouptut;
            auto memory_time = execute_command(folder + "evaluation/tmp/" + "query_output.txt", command);
            memories.push_back( std::get<0>(memory_time));
            times.push_back( std::get<1>(memory_time));

        }
        assert(times.size() == number_of_queries);

        double memory = std::accumulate(memories.begin(), memories.end(), 0.0) / memories.size();
        mean = memory;
        double std_memory = std::sqrt(std::accumulate(memories.begin(), memories.end(), 0.0, variance_func));

        double time =std::accumulate(times.begin(), times.end(), 0.0) / times.size();
        mean = time;
        double std_time = std::sqrt(std::accumulate(times.begin(), times.end(), 0.0, variance_func));
        std::cout << "std time" << std_time << std::endl;

        return std::make_tuple(memory, time);
    }

    // future TODO check if the search results are correct using a test function
}








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

long long convert_to_numeric(std::string size_with_suffix) {
    char suffix = size_with_suffix.back();
    size_with_suffix.pop_back();
    long long size;
    std::istringstream(size_with_suffix) >> size;

    try {
        switch (suffix) {
            case 'K':
                return size * 1000;
            case 'M':
                return size * 1000000;
            case 'G':
                // Check if comma is present before replacing
                if (size_with_suffix.find(',') != std::string::npos) {
                    size_with_suffix.replace(size_with_suffix.find(','), 1, ".");
                }
                return stod(size_with_suffix) * 1000000000;
            default:
                return size;
        }
    }
    catch (const std::exception& e) {
        // Handle the exception (print an error message, etc.)
        std::cerr << "Error occurred: " << e.what() << std::endl;
        return 0; // Return a default value or handle the error appropriately
    }
}

int file_size(std::string filename_index){
            // measure the size of the uncompressed index.
        std::ifstream file_index(filename_index, std::ifstream::binary);
        file_index.seekg(0, file_index.end);
        int fileSize = file_index.tellg(); // number of bytes.
        file_index.seekg(0, file_index.beg);
        std::cout << "filesize CPP: " << fileSize << std::endl;
        return fileSize;
}


int file_size_bash(std::string filename_index){
            // measure the size of the uncompressed index.
            system(("du -s " + filename_index + " | cut -f1 > temp_file.txt" ).c_str()); // returns KB
            std::string result;
            std::ifstream temp_file("temp_file.txt");
            if (temp_file) {
                std::getline(temp_file, result);
                temp_file.close();
            }

            // Remove newline characters from the end of the result string
            if (!result.empty() && result.back() == '\n') {
                result.pop_back();
            }

            // Remove the temporary file
            std::remove("temp_file.txt");


        std::cout << "filesize system: " << result << std::endl;
        // TODO use also the second method to obtain the filesize.
        //std::cout << "filesize system: "  << system(("du -sh " + filename_index + " | cut -f1)").c_str()) << std::endl;
        //auto fileSize = convert_to_numeric(result);
        int fileSize = std::stod(result);
        return fileSize;
}


template <typename T>
std::string outstring(std::vector<T> out_array){
    std::string out = "";
    for (auto &item: out_array) { std::string x = std::to_string(item); out = out + x + ", "; }
    return out;
}

// SAVE RESULTS
// Saving results to a python library, that can be read in to make figures.
int write_to_python(std::string python_filename,
                    std::vector<double> time_insertion, std::vector<double> time_query,
                    std::vector<double> memory_insertion, std::vector<double> memory_query,
                    std::vector<int> size_index, std::vector<int> size_index_bash, std::vector<int> rebuild){
    std::ofstream o(python_filename + ".py"); // "here.py");//
    o << "time_insertion = [" << outstring(time_insertion) << "]" <<std::endl;
    o << "time_query = [" << outstring(time_query) << "]" <<std::endl;
    o << "memory_insertion = [" << outstring(memory_insertion) << "]" <<std::endl;
    o << "memory_query = [" << outstring(memory_query) << "]" <<std::endl;
    o << "size_index = [" << outstring(size_index) << "]" <<std::endl;
    o << "size_index_bash = [" << outstring(size_index) << "]" <<std::endl;
    o << "rebuild = [" << outstring(rebuild) << "]" <<std::endl;

    o.close();
    return 0;
}



//////////////////////////////////////////////////

int main_insert_ub(){
    std::string folder = std::filesystem::current_path(); folder += "/";
    std::cout << "Give the test_folder, should be located within the evaluation folder: ";
    std::string input_test; std::cin >> input_test;
    std::cout << "Give the tmp folder, should be located within the evaluation folder: ";
    std::string tmp_folder; std::cin >> tmp_folder;
// PARAMETERS
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";"evaluation.index"; // this could best be an index without empty bins.

    std::string filename_executable = folder +  "raptor";

    // output
    std::string python_filename = folder + "evaluation/" + input_test + "/results/";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";
    double sample_percentage = 0.001;

    //std::filesystem::remove_all("tmp");
    system("mkdir evaluation");
    system(("mkdir " + python_filename).c_str());
    std::string filename_queries_existing = folder + "evaluation/" + tmp_folder + "/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/" + tmp_folder + "/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/" + tmp_folder).c_str());
    std::string existing_filenames_building = folder + "evaluation/" + tmp_folder + "/" + "existing_filenames.txt";
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
    std::cout << python_filename <<std::flush;



    auto result = create_query_file(existing_filenames_building, filename_queries_existing, sample_percentage, folder, tmp_folder);
    std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);



for (auto insertion_method: {"find_ibf_size_splitting", "find_ibf_idx_traverse_by_fpr", "naive", "find_ibf_idx_traverse_by_similarity",  "find_ibf_idx_ibf_size",   }){
    std::string filename_index =  folder + "evaluation/" + tmp_folder + "/" +insertion_method+ "_" +  filename_index_original;
    int number_of_files = extract_filenames(existing_filenames_building).size(); // TODO write a function get number of lines instead.
    std::cout << std::endl << std::endl  << insertion_method << std::endl  << std::flush;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file
    // TODO why is the query time so much higher for the FPR method?

    // Result vectors
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds;
    // check time and size before updating
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));
    size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));
    int counter=0;
    for (const std::string& user_bin_filename : user_bin_filenames) {
        // TODO change such that we can include a batch size.
        counter += 1;
        system(("mkdir " + tmp_folder + "").c_str()); //
        number_of_files += 1;
        std::cout << "Method: " << insertion_method << ", with current bin number: " << number_of_files << std::endl << std::flush;
        // store filename to a temporary file.
        size_t lastSlashPos = user_bin_filename.find_last_of("/");
        std::string lastPart = user_bin_filename.substr(lastSlashPos + 1);
        std::string tmp_filename = folder + "evaluation/" + tmp_folder + "/" + lastPart;
        std::ofstream tmp_insertion_filename(tmp_filename);
        tmp_insertion_filename << user_bin_filename;
        tmp_insertion_filename.close();
        // update all paths and queries
        if (generate_reads){
                auto result = create_query_file(existing_filenames_building, filename_queries_existing, sample_percentage, folder);
        }else{
                    write_to_fasta(filename_queries, user_bin_filename, sample_percentage); // all queries together
        }
        // because otherwise you might measure the time difference by parallelization
        // or increase query size and put threads to 0, but that does not make it very comparable to mantis.

        std::string tmp_query_filename = folder + "evaluation/" + tmp_folder + "/" + lastPart;
        tmp_query_filenames.push_back(tmp_query_filename);
        //existing_filenames.push_back(user_bin_filename);

        std::tuple<int, double, bool, bool> memory_time_insertion;
        if (insertion_method != "naive"){
            memory_time_insertion = insert_ub(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, insertion_method, tmp_folder); //  measure insertion time and memory
        }else{
            write_to_txt(existing_filenames_building, user_bin_filename);
            memory_time_insertion = rebuild_index(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, existing_filenames_building, tmp_folder); //  measure insertion time and memory

        }
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        if (std::get<3>(memory_time_insertion) and not std::get<2>(memory_time_insertion)){ // better distinghuis between partial and full rebuild.  and not std::get<2>(memory_time_insertion)
            rebuilds.push_back(2); // 2 inidicates a partial rebuild.
        }else{
            rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        }

        // when not using "all_at_once", then empty tmp_query_filename
        // write_to_fasta(tmp_query_filename, user_bin_filename, sample_percentage); // with a single query

        auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));

        size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));
        if (counter%50==0){
            std::cout << "saving intermediate results" <<std::endl;
                    write_to_python(python_filename + insertion_method, time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
        }

    }

    write_to_python(python_filename + insertion_method, time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
}
// TODO test if all files have been added. use load_hibf
// TODO built in a better way to ensure that each insertion/query was without errors
// TODO find index by size doesn't seem to go so well.
// delete tmp files.
    return 0;
}

//INSERT SEQUENCES
//1. Insert single bin and measure insertion time.
std::tuple<int, double, bool, bool>  insert_sequences(std::string filename_ub, std::string filename_index,
                                   std::string filename_executable, std::string sketch_directory,
                                   std::string folder){
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
    std::string ouptut_file = folder + "evaluation/tmp/" + "insertion_output.txt";
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
    std::cout << "Give the test_folder, should be located within the evaluation folder: ";
    std::string input_test; std::cin >> input_test;
// PARAMETERS
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";"evaluation.index"; // this could best be an index without empty bins.

    std::string filename_executable = folder +  "raptor";

    // output
    system(("mkdir " + folder + "evaluation/" + input_test + "/results/").c_str());
    std::string python_filename = folder + "evaluation/" + input_test + "/results/results";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";
    double sample_percentage = 0.001;

    //std::filesystem::remove_all("tmp");
    system("mkdir evaluation");
    system("mkdir evaluation/results");
    std::string filename_queries_existing = folder + "evaluation/tmp/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/tmp/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/tmp").c_str());
    std::string existing_filenames_building = folder + "evaluation/tmp/" + "existing_filenames.txt";
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
    std::cout << python_filename <<std::flush;



    auto result = create_query_file(existing_filenames_building, filename_queries_existing, sample_percentage, folder);
    int number_of_files = extract_filenames(existing_filenames_building).size(); std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::cout << "number of lines: " << number_of_files << std::endl;
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);

    std::string filename_index =  folder + "evaluation/tmp/" +"sequence_insertions"+ "_" +  filename_index_original;

    std::cout << std::endl << std::endl  << "sequence_insertions" << std::endl  << std::flush;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file

    // Result vectors
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds;
    // check time and size before updating
    auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames); //  measure query times
    memory_query.push_back(std::get<0>(memory_time_queries));
    time_query.push_back(std::get<1>(memory_time_queries));
    size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));

        system("mkdir tmp"); //


        //create file with sequence content to insert "bin_0000_insertsequences"
        // pick next 10 lines from bin_0000

    std::string file_to_insert_from = "example_data/64/bins/bin_bonus.fasta";
    std::ifstream inputFile(file_to_insert_from);  // Open input file
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open input file." << std::endl;
        return 1;
    }
    std::string insert_to_ub = "example_data/1024/bins/bin_0000.fasta";
    //  alternatively randomly sample from the existing user bins.
//        std::random_device rd;
//        std::mt19937 rng(rd());
//        std::uniform_int_distribution<size_t> dist(0, user_bin_filenames.size() - 1);
//        size_t random_index = dist(rng);
//        std::string insert_to_ub = user_bin_filenames[random_index];
    std::string insert_to_ub_paths = "evaluation/tmp/insert_to_ub_paths.txt";
    std::ofstream pathFile(insert_to_ub_paths);
    pathFile << insert_to_ub << std::endl;
    pathFile.close();
    std::string filename_new_sequences = insert_to_ub.substr(0, insert_to_ub.find_last_of('.')) + "_insertsequences"
            + insert_to_ub.substr(insert_to_ub.find_last_of('.'));

    int fileCount = 0;
    int lineCount = 0;
    std::string line;

    std::ofstream outputFile(filename_new_sequences, std::ios_base::trunc);  // Open output file in append mode
    outputFile << ">Header1\n";
    outputFile.close();
    std::system(("cp " + insert_to_ub + " " +  insert_to_ub + "_original").c_str() ); // if only 1 file is updated.

    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '>') {
            continue;  // Skip empty lines or lines starting with '>'
        }

        if (lineCount % 300 == 0 and lineCount!= 0) { // 10 lines will be inserted a time.
        std::system(("cat " + filename_new_sequences + " >> " + insert_to_ub ).c_str());


         std::tuple<int, double, bool, bool> memory_time_insertion = insert_sequences(insert_to_ub_paths, filename_index,
                               filename_executable, sketch_directory,
                               folder);
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        if (std::get<3>(memory_time_insertion) and not std::get<2>(memory_time_insertion)){
            rebuilds.push_back(2); // 2 inidicates a partial rebuild.
        }else{
            rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        }
        auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames); //  measure query times
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));

        size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));

        std::remove(filename_new_sequences.c_str());         // empty the file outputFileName, such that 10 new sequences can be inserted.
        std::ofstream outputFile(filename_new_sequences);  // Open output file in append mode
        outputFile << ">Header\n";     outputFile.close();

        }

        if (lineCount%500==0){
            std::cout << "saving intermediate results" <<std::endl;
                    write_to_python(python_filename + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
        }

        std::ofstream outputFile(filename_new_sequences, std::ios_base::app);  // Open output file in append mode
        outputFile << line << std::endl; outputFile.close();
        lineCount++;
    }

    inputFile.close();

    write_to_python(python_filename + "_insertsequences", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
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

// PARAMETERS
    std::string insertion_paths = folder + "evaluation/" + input_test + "/insertion_paths.txt"; //"update_bin_paths_multiple.txt";
    std::string all_paths = folder + "evaluation/" + input_test + "/existing_paths.txt";;//"half_of_bin_paths.txt"; //"all_bin_paths.txt";// existing bin paths, used for querying all bins.
    std::string filename_index_original = "hibf.index";"evaluation.index"; // this could best be an index without empty bins. // start with an index or around 100 bins.
    std::string filename_executable = folder +  "raptor";

    // output
    std::string python_filename = folder + "evaluation/" + input_test + "/results/";
    std::string sketch_directory = folder +  "chopper_sketch_sketches";
    double sample_percentage = 0.001;

    //std::filesystem::remove_all("" + tmp_folder + "");
    system("mkdir evaluation");
    system("mkdir evaluation/results");
    std::string filename_queries_existing = folder + "evaluation/" + tmp_folder + "/" + "queries_original.fasta";
    std::string filename_queries = folder + "evaluation/" + tmp_folder + "/" + "queries.fasta";
    system(("mkdir " +folder + "evaluation/" + tmp_folder + "").c_str());
    std::string existing_filenames_building = folder + "evaluation/" + tmp_folder + "/" + "existing_filenames.txt";
    std::remove(existing_filenames_building.c_str()); // delete file first
    system(("yes | cp -f " + all_paths + " " + existing_filenames_building).c_str()); //make a copy of the file
    //todo add endl to existing filenames.

    ////////////
    system("cp ../_deps/raptor_chopper_project-src/build/bin/chopper chopper");

    std::ofstream config("evaluation/results/configuration");
    config  << "folder: " << folder <<std::endl;
    config << "tmp_folder: " << tmp_folder <<std::endl;
    config << "kmer_size: " << kmer_size <<std::endl;
    config << "window_size: " << window_size <<std::endl;
    config << "fpr: " << fpr <<std::endl;
    config << "num_hash_functions: " << num_hash_functions <<std::endl;
    config << "kmer_size: " << kmer_size <<std::endl;
    config << "query_threads: " << query_threads <<std::endl;
    config << "query_errors: " << query_errors <<std::endl;
    config << "generate_reads: " << kmer_size <<std::endl;
    config << "number_of_reads: " << kmer_size <<std::endl;
    config << "read_length: " << kmer_size <<std::endl<<std::endl;

    config << "insertion_paths: " << insertion_paths <<std::endl;
    config << "all_paths " << all_paths <<std::endl;
    config << "filename_index_original: "<< filename_index_original <<std::endl;
//    config << filename_index <<std::endl;
    config << "existing_filenames_building:" << existing_filenames_building <<std::endl;
    config << "filename_queries_existing: " << filename_queries_existing <<std::endl;
    config << "filename_queries: " << filename_queries <<std::endl;
    config << "sketch_directory: " << sketch_directory <<std::endl;
    config << "python_filename: " << python_filename <<std::endl;
    config.close();


    auto result = create_query_file(existing_filenames_building, filename_queries_existing, sample_percentage, folder, tmp_folder);
    system(("yes | cp -rf " + filename_queries_existing + " " + filename_queries).c_str()); //make a copy of the file

    int number_of_files = extract_filenames(existing_filenames_building).size(); std::vector<std::string> tmp_query_filenames = std::get<1>(result);
    std::vector<std::string> user_bin_filenames = extract_filenames(insertion_paths);
//o	Iteratively delete user bins?
//o	Misschien in combinatie met seq insertions; insert delete insert etc. -> laten zien dat het niet groeit.
// UB insertion
// UB deletion ; all of the same size.
    std::string tmp_filename = folder + "evaluation/" + tmp_folder + "/tmp_input_file.txt";
    std::vector<double> time_insertion, time_query, memory_insertion, memory_query;
    std::vector<int> size_index, size_index_bash, rebuilds, del_or_insert;
    std::string insertion_paths_copy = "evaluation/" + tmp_folder + "/insertion_paths_copy.txt ";
    system(("yes | cp -f " + insertion_paths + " " + insertion_paths_copy).c_str()); //make a copy of the file

    std::string filename_index =  folder + "evaluation/" + tmp_folder + "/" "deletions" + "_" +  filename_index_original;
    system(("yes | cp -rf " + folder + "evaluation/" + input_test + "/" + filename_index_original + " " + filename_index).c_str()); //make a copy of the file


 // delete insert delete delete insert insert
 int user_bin_filenames_counter =0;
 for (int number_of_operatons = 0; number_of_operatons < 1000; number_of_operatons++){
     if (user_bin_filenames_counter + number_of_operatons < user_bin_filenames.size()){
     for (int _ = 0; _ < number_of_operatons; _++){
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
            auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames, true, tmp_folder); //  measure query times
            memory_query.push_back(std::get<0>(memory_time_queries));
            time_query.push_back(std::get<1>(memory_time_queries));
            size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));
         }


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

        auto memory_time_insertion = insert_ub(tmp_filename, filename_index,  filename_executable, sketch_directory, folder, "find_ibf_idx_traverse_by_fpr", tmp_folder); //  measure insertion time and memory
        memory_insertion.push_back(std::get<0>(memory_time_insertion));
        time_insertion.push_back(std::get<1>(memory_time_insertion));
        rebuilds.push_back((int) std::get<2>(memory_time_insertion));
        del_or_insert.push_back(0);
        auto memory_time_queries = query_all_ubs(filename_queries, filename_index, filename_executable, number_of_files, folder, tmp_query_filenames); //  measure query times
        memory_query.push_back(std::get<0>(memory_time_queries));
        time_query.push_back(std::get<1>(memory_time_queries));
        size_index.push_back(file_size(filename_index));
    size_index_bash.push_back(file_size_bash(filename_index));

}
     }
     // expected result: rebuild after some time because merged bins grow in FPR.
     // before that, querying might take longer because of false hits. (depends on how much you query and if you insert the UBs that you deleted, how similart they are)
     write_to_python(python_filename + "deletions", time_insertion, time_query, memory_insertion, memory_query, size_index, size_index_bash, rebuilds); //  write result vectors to python file.
//TODO also write del_or_inserts to python file.
 }
}


int main(){
    std::cout<< "seq/ub/del" ;
    std::string test_type; std::cin >> test_type;
    if (test_type == "ub"){
        main_insert_ub();
    }else if (test_type == "seq"){
        main_insert_seq();
    }else if (test_type == "del"){
        main_del_seq();
    }
    return 0;
}