#include <iostream>
#include <fstream>
#include <string>

int main() {
    // Open input files
    std::ifstream file1("./rv_left.dat");
    std::ifstream file2("./rv_right.dat");
    
    // Check if files are open
    if (!file1.is_open()) {
        std::cerr << "Error: Could not open rv_left.dat" << std::endl;
        return EXIT_FAILURE;
    }
    
    if (!file2.is_open()) {
        std::cerr << "Error: Could not open rv_right.dat" << std::endl;
        return EXIT_FAILURE;
    }

    // Open output file
    std::ofstream combined_file("rv_start.dat");
    
    if (!combined_file.is_open()) {
        std::cerr << "Error: Could not create rv_start.dat" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Write contents of file1 to combined_file
    combined_file << file1.rdbuf();
    
    // Write contents of file2 to combined_file
    combined_file << file2.rdbuf();
    
    // Close the files
    file1.close();
    file2.close();
    combined_file.close();
    
    return EXIT_SUCCESS;
}
