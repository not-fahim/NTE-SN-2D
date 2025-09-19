#include "geo_mat_parser.hpp"
#include <iostream>
#include <iomanip> // For formatting the output

// A helper function to print the results
void print_results(const input_class& data) {
    if (data.name == "nan") {
        std::cout << "Failed to read input file." << std::endl;
        return;
    }

    std::cout << "--- Successfully Parsed Data ---" << std::endl;
    std::cout << "Name: " << data.name << std::endl;
    std::cout << "SN Value: " << data.S_n << std::endl;

    std::cout << "\nCellX Data:" << std::endl;
    for (const auto& cell : data.cellX) {
        std::cout << "  Dimension: " << cell.cell_dim << ", Mesh Number: " << cell.mesh_num << std::endl;
    }

    std::cout << "\nCellY Data:" << std::endl;
    for (const auto& cell : data.cellY) {
        std::cout << "  Dimension: " << cell.cell_dim << ", Mesh Number: " << cell.mesh_num << std::endl;
    }

    std::cout << "\nBoundary Conditions (BC): ";
    for (double val : data.BC) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    std::cout << "\nMaterial Map:" << std::endl;
    for (const auto& row : data.matmap_cell) {
        std::cout << "  ";
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------------------" << std::endl;
}


int main() {
    // Call your function to parse the file
    input_class parsed_data = read_input_file();

    // Print the contents of the returned struct
    print_results(parsed_data);
    
    return 0;
}