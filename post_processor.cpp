#include "post_processor.hpp"
post_processor_class::post_processor_class(geometry_class const &geometry_object,
                         std::vector<std::vector<std::vector<double>>> const &flux_input, 
                         input_class const &input_object)
{

    int mesh_num_x = geometry_object.Nx;
    int mesh_num_y = geometry_object.Ny;
    int refinement = input_object.refinement;
    int groups = geometry_object.groups;

    int cell_num_x = input_object.cellX_vector.size();
    int cell_num_y = input_object.cellY_vector.size();
    int sub_cell_x_num, sub_cell_y_num;

    int iii=0, jjj=0; // keeps track of sub-cell indices in the universal mesh
    flux_g = std::vector<std::vector<std::vector<double>>> (groups, std::vector<std::vector<double>>(cell_num_x, std::vector<double>(cell_num_y,0.0)));
    fission_density_g = std::vector<std::vector<std::vector<double>>> (groups, std::vector<std::vector<double>>(cell_num_x, std::vector<double>(cell_num_y,0.0)));
    for(int g=0; g<groups; g++)
    {
        for(int i=0; i<cell_num_x; i++)
        {
            sub_cell_x_num = input_object.cellX_vector[i].mesh_num * refinement;
            for(int ii=0; ii<sub_cell_x_num; ii++)
            {
                for(int j=0; j<cell_num_y; j++)
                {   
                    jjj = 0;  
                    sub_cell_y_num = input_object.cellY_vector[j].mesh_num * refinement;
                    for(int jj=0; jj<sub_cell_y_num; jj++)
                    {
                        flux_g[g][i][j] += flux_input[g][iii][jjj];
                        fission_density_g[g][i][j] += flux_input[g][iii][jjj]*geometry_object.nu_sigmaf(iii, jjj, g)   ;
                        jjj++;
                    }
                }
                iii++;
            }
        }
    }
}