#include "output_handling.hpp"

/**
 * @brief Writes the final scalar flux and geometry information to an output file.
 */
void write_flux_output(
    const input_class& input,
    const geometry_class& geometry,
    const std::vector<std::vector<std::vector<double>>>& flux)
{
    // Create the output filename
    std::string filename = "result_" + input.name + ".flux";
    std::ofstream outfile(filename);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    // Set precision for floating point numbers
    outfile << std::fixed << std::setprecision(8);

    // Write the header comment
    outfile << "% Flux array organization: [group][i_x][j_y]" << std::endl;
    outfile << "% Flux matrices below are printed with i_x (columns) and j_y (rows)." << std::endl;
    outfile << std::endl;

    // Write cellx array (x_i midpoints)
    outfile << "cellx";
    for (int i = 0; i < geometry.Nx; ++i) {
        outfile << " " << geometry.x_i[i];
    }
    outfile << std::endl;

    // Write celly array (y_j midpoints)
    outfile << "celly";
    for (int j = 0; j < geometry.Ny; ++j) {
        outfile << " " << geometry.y_j[j];
    }
    outfile << std::endl;

    // Write refinement
    outfile << "refinement " << input.refinement << std::endl;
    outfile << std::endl; // Add a newline for spacing

    // Write the flux for each group
    for (int g = 0; g < geometry.groups; ++g) {
        outfile << "phi" << (g + 1) << std::endl;

        // Print the 2D flux matrix for group g
        // We print [j][i] to get Ny rows and Nx columns
        for (int j = 0; j < geometry.Ny; ++j) {
            for (int i = 0; i < geometry.Nx; ++i) {
                // Accessing flux as [g][i][j]
                outfile << flux[g][i][j] << ((i == geometry.Nx - 1) ? "" : " ");
            }
            outfile << std::endl;
        }
        outfile << std::endl; // Add a newline for spacing between groups
    }

    outfile.close();
    std::cout << "Successfully wrote flux output to " << filename << std::endl;
}
