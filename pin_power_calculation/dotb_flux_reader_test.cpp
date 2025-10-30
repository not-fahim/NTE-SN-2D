/**
 * @file main.cpp
 * @brief Test program for the read_dotb_flux_ouput function.
 * * This program reads a .flux file specified by 'filename' and prints the
 * contents of the resulting doutb_flux_output_class struct to the console
 * for verification.
 * * @note You must have the file "result_34x34-Bare-3x3samepinxs.flux" in the
 * same directory as the compiled executable for this test to work.
 * * @compile
 * To compile this test program, you need to link both this file and your
 * implementation file (pin_power.cpp) together.
 * * Example (using g++):
 * g++ main.cpp pin_power.cpp -o test_parser -std=c++11
 * * @run
 * After compiling, run the executable:
 * ./test_parser
 */

#include <iostream>
#include <string>
#include <vector>
#include <iomanip> // For std::setw, std::fixed, std::setprecision

// Include the header file that defines doutb_flux_output_class
// and declares the read_dotb_flux_ouput function.
#include "pin_power.hpp" 

/**
 * @brief A helper template function to print the contents of a vector.
 * @tparam T The type of elements in the vector.
 * @param name A descriptive name for the vector (e.g., "pin_num").
 * @param vec The vector to be printed.
 */

void print3DVector(const std::vector<std::vector<std::vector<double>>>& vec) {
    // Loop through the first dimension (index i)
    for (size_t i = 0; i < vec.size(); ++i) {
        // Print the header for this 2D "slice"
        std::cout << "vector[" << i << "]:" << std::endl;

        // Loop through the second dimension (index j)
        for (size_t j = 0; j < 10; ++j) {
            
            // Loop through the third dimension (index k)
            // This loop prints all elements of vec[i][j] on one line
            for (size_t k = 0; k < 10; ++k) {
                std::cout << vec[i][j][k] << " ";
            }
            // Move to the next line after printing a full 1D row
            std::cout << std::endl; 
        }
        // Optional: Add an extra newline between 2D slices for clarity
        // cout << endl; 
    }
}

void print2DVector(const std::vector<std::vector<double>> &vec) {

        for (size_t j = 0; j < 10; ++j) {
            
            // Loop through the third dimension (index k)
            // This loop prints all elements of vec[i][j] on one line
            for (size_t k = 0; k < 10; ++k) {
                std::cout << vec[j][k] << " ";
            }
            // Move to the next line after printing a full 1D row
            std::cout << std::endl; 
        }

}

template<typename T>
void print_vector(const std::string& name, const std::vector<T>& vec) {
    std::cout << "  " << std::setw(12) << std::left << (name + ":") << "[ ";
    // Print only the first 10 elements for brevity if the vector is large
    size_t limit = std::min(vec.size(), (size_t)10);
    for (size_t i = 0; i < limit; ++i) {
        std::cout << vec[i] << " ";
    }
    if (vec.size() > 10) {
        std::cout << "... (" << vec.size() << " total elements)";
    }
    std::cout << "]" << std::endl;
}

int main() {
    // Define the filename of the test file.
    // This file must be in the same directory as the executable.
    std::string filename = "result_34x34-Bare-3x3samepinxs.flux";

    std::cout << "Attempting to read file: " << filename << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    // Call the function to test
    doutb_flux_output_class flux_data = read_dotb_flux_ouput(filename);

    // Check for the error flag (name == "nan")
    if (flux_data.name == "nan") {
        std::cerr << "Test FAILED: Function returned error flag." << std::endl;
        std::cerr << "Ensure the file '" << filename << "' exists and has read permissions." << std::endl;
        return 1; // Exit with an error code
    }

    // If no error, print all the loaded data
    std::cout << "Test PASSED. File read successfully." << std::endl;
    std::cout << "Data loaded into struct:" << std::endl;
    
    std::cout << "  " << std::setw(12) << std::left << "name:" << flux_data.name << std::endl;
    std::cout << "  " << std::setw(12) << std::left << "groups:" << flux_data.groups << std::endl;
    std::cout << "  " << std::setw(12) << std::left << "pin_calc:" << (flux_data.pin_calc ? "true" : "false") << std::endl;
    std::cout << "  " << std::setw(12) << std::left << "pin_pitch:" << std::fixed << std::setprecision(8) << flux_data.pin_pitch << std::endl;

    // Print vector data using the helper function
    print_vector("pin_num", flux_data.pin_num);
    print_vector("mesh_xi", flux_data.mesh_xi);
    print_vector("mesh_yj", flux_data.mesh_yj);
    print_vector("mesh_delxi", flux_data.mesh_delxi);
    print_vector("mesh_delyj", flux_data.mesh_delyj);

    // Print the 3D flux vector (mesh_flux[group][y][x])
    std::cout << "\n  --- Mesh Flux Data (mesh_flux) ---" << std::endl;
    std::cout << "  (Showing size and first 5x5 elements of each group)" << std::endl;
    print3DVector(flux_data.mesh_flux);
    std::cout << "\n  --- Mesh Fission Density Data (mesh_fission_density) ---" << std::endl;
    print2DVector(flux_data.mesh_fission_density);
    

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Test finished." << std::endl;

    return 0; // Exit successfully
}
