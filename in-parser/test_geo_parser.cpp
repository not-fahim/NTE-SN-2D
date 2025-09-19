#include "geo_mat_parser.hpp"
#include <iostream>
#include <iomanip> // For formatting the output

// Updated function to print all parsed data, including materials
void print_results(const input_class& data) {
    if (data.name == "nan") {
        std::cout << "Failed to read input file." << std::endl;
        return;
    }
    if (data.name == "no_mat") {
        std::cout << "after reading geometry, Failed to read material input file." << std::endl;
        return;
    }

    std::cout << "--- Successfully Parsed Geometry Data ---" << std::endl;
    std::cout << "Name: " << data.name << std::endl;
    std::cout << "Material File: " << data.xsfile << std::endl;
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
    std::cout << "---------------------------------------" << std::endl;

    // --- New section to print material data ---
    std::cout << "\n--- Successfully Parsed Material Data ---" << std::endl;
    std::cout << std::fixed << std::setprecision(8); // For better float output

    for (const auto& mat : data.mat_vector) {
        std::cout << "\n## Material: " << mat.name << " ##" << std::endl;

        std::cout << "Total (sigma_t):      ";
        for (double val : mat.sigma_t) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Nu-Fission (nu_sigma_f):";
        for (double val : mat.nu_sigma_f) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Fission Spectrum (chi): ";
        for (double val : mat.chi) std::cout << val << "  ";
        std::cout << std::endl;

        std::cout << "Scattering Matrix (sigma_s):" << std::endl;
        for (const auto& row : mat.sigma_s) {
            std::cout << "  ";
            for (double val : row) {
                std::cout << val << "  ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "---------------------------------------" << std::endl;
}

int main() {
    // Call your function to parse the file
    input_class parsed_data = read_input_file();

    // Print the contents of the returned struct
    print_results(parsed_data);
    
    return 0;
}