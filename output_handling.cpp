#include "output_handling.hpp"

void print_banner(std::ofstream &outputfile) {
    // Copy and paste the text from one of the banners below into this string
    std::string banner_text = R"(
%********************************************************
%*                 .-~*~~~*~-.                          *
%*                .-~~~~~~~~~-.                         *
%*               /  X       X  \                        *
%*              |    .-----.    |                       *
%*               \  '_______'  /                        *
%*                `-.........-'                         *
%* *          ____   ___  ____ ____                     *
%*           |  _ \ / _ \_   _| __ )                    *
%*           | | | | | | || | |  _ \                    *
%*           | |_| | |_| || | | |_) |                   *
%*           |____/ \___/ |_| |____/                    *
%********************************************************
%*      Discrete Ordinate with Two Braincells           *
%********************************************************
    )";
    outputfile << banner_text << std::endl;
};


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

    if (!outfile.is_open()) 
    {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return;
    }

    // Set precision for floating point numbers
    outfile << std::fixed << std::setprecision(8);
    print_banner(outfile);
    // Write the header comment
    outfile<<"%%%%Discrete Ordinates with Two BrainCells%%%%"<< std::endl;
    outfile<<"%%%%%%%%%Pin Power Calculation Output%%%%%%%%%"<<std::endl << std::endl;
    outfile << "% Flux array organization: [group][i_x][j_y]" << std::endl;
    outfile << "% Flux matrices below are printed with i_x (columns) and j_y (rows)." << std::endl;
    outfile << "% an input parser reads the first row, stores in a vector and so on. so vector of vectors will have access index [i_x][i_y]" << std::endl;

    outfile << std::endl;

    if(input.pin_calc == false) 
    {
        outfile << "% this output was not generated for pin power calculation from DOTB"<<std::endl;
        outfile<<"pincalc false"<<std::endl;
        outfile << std::endl;
        outfile << "%% if you want to calculate pin flux, modify orginal dotb input problem file and rerun the solver" << std::endl;
        outfile << std::endl;
    } 
    else
    {
        outfile << "% this output was generated for pin power calculation from DOTB"<<std::endl;
        outfile << "pincalc true" << std::endl;
        outfile<< "pitch "<<input.pin_pitch << " cm" << std::endl;
        outfile << "pins " << input.pin_num[0] <<"\t" << input.pin_num[1] << std::endl;
    }
    outfile << std::endl;
    // Write cellx array (x_i midpoints)
    outfile<< "% the cell centers used in the mesh grid in solver"<<std::endl;
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
    if(input.pin_calc == true)
    {
        outfile << std::endl<<std::endl;
        
        outfile << "% array of length of each cell in mesh used in the solver & is useful for pin power calc"<<std::endl;
        //write delx
        outfile << "delx";
        for (int i = 0; i < geometry.Nx; ++i) {
            outfile << " " << geometry.delx_i[i];
        }
        outfile << std::endl;
        //write dely
        outfile << "dely";
        for (int j = 0; j < geometry.Ny; ++j) {
            outfile << " " << geometry.dely_j[j];
        }
    }
    outfile << std::endl;
    outfile << std::endl; // Add a newline for spacing

    // Write the flux for each group
    for (int g = 0; g < geometry.groups; ++g) 
    {
        outfile << "phi" << (g + 1) << std::endl;

        // Print the 2D flux matrix for group g
        // We print [j][i] to get Ny rows and Nx columns
        for (int j = 0; j < geometry.Ny; ++j) 
        {
            for (int i = 0; i < geometry.Nx; ++i) 
            {
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
