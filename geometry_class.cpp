#include "geometry_class.hpp"

std::vector<double> range(double min, double max, size_t N)
{
    //returns a vector of N equally spaced doubles including min  max
    std::vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}

std::vector<double> cell_midpoints(std::vector<double> x_i_minus_half, std::vector<double> x_i_plus_half)
{
    std::vector<double> cell_midpoints;
    for(int i=0; i<size(x_i_minus_half); i++)
    {
        cell_midpoints.push_back((x_i_minus_half[i] + x_i_plus_half[i])/2);
    }
    return cell_midpoints;
}

geometry_class::geometry_class(input_class const &input_object)
{
    //vectors to work inside function. thus have "in" in the name
    std::vector<double> cell_xi_minus_half, cell_yj_minus_half, cell_xi_plus_half, cell_yj_plus_half;
    std::vector<double> cell_xi, cell_yj;
    std::vector<double> cell_i, cell_j;
    int cell_Nx, cell_Ny;
    double cell_hx, cell_hy, cellx_dim, celly_dim;

    //materials and groups
    materials = input_object.mat_vector;
    groups = materials[0].sigma_t.size();
    
    //boundary condition
    alpha_B = input_object.BC[0];
    alpha_T = input_object.BC[1];
    alpha_L = input_object.BC[2];
    alpha_R = input_object.BC[3];

    //meshing
    double x_in, x_out, y_in, y_out;

    for(int cellx = 0; cellx < input_object.cellX_vector.size(); cellx++)
    {
        if(cellx == 0)
        {
            x_in = 0;
        }
        else
        {
            x_in = x_out;
        }
        cellx_dim = input_object.cellX_vector[cellx].cell_dim;
        x_out = x_in + cellx_dim;

        cell_Nx = input_object.refinement*input_object.cellX_vector[cellx].mesh_num;
        cell_hx = cellx_dim/cell_Nx;

        cell_xi_minus_half = range(x_in, x_out-cell_hx, cell_Nx);
        cell_xi_plus_half = range(x_in+cell_hx, x_out, cell_Nx);
        cell_xi = cell_midpoints(cell_xi_minus_half, cell_xi_plus_half);

        for(int i= 0; i< cell_xi.size(); i++)
            cell_i.push_back(cellx);

        x_i.insert(x_i.end(), cell_xi.begin(), cell_xi.end());
        x_i_minus_half.insert(x_i_minus_half.end(), cell_xi_minus_half.begin(), cell_xi_minus_half.end());
        x_i_plus_half.insert(x_i_plus_half.end(), cell_xi_plus_half.begin(), cell_xi_plus_half.end());


    }

    for(int celly = 0; celly < input_object.cellY_vector.size(); celly++)
    {
        if(celly == 0)
        {
            y_in = 0;
        }
        else
        {
            y_in = y_out;
        }
        celly_dim = input_object.cellY_vector[celly].cell_dim;
        y_out = y_in + celly_dim;

        cell_Ny = input_object.cellY_vector[celly].mesh_num * input_object.refinement;
        cell_hy = celly_dim/cell_Ny;

        cell_yj_minus_half = range(y_in, y_out-cell_hy, cell_Ny);
        cell_yj_plus_half = range(y_in+cell_hy, y_out, cell_Ny);
        cell_yj = cell_midpoints(cell_yj_minus_half, cell_yj_plus_half);

        for(int j= 0; j< cell_yj.size(); j++)
            cell_j.push_back(celly);

        y_j.insert(y_j.end(), cell_yj.begin(), cell_yj.end());
        y_j_minus_half.insert(y_j_minus_half.end(), cell_yj_minus_half.begin(), cell_yj_minus_half.end());
        y_j_plus_half.insert(y_j_plus_half.end(), cell_yj_plus_half.begin(), cell_yj_plus_half.end());

    }

    for (int i = 0; i < x_i.size(); i++)
    {
        delx_i.push_back(x_i_plus_half[i] - x_i_minus_half[i]);
    }

    for (int j = 0; j < y_j.size(); j++)
    {
        dely_j.push_back(y_j_plus_half[j]- y_j_minus_half[j] );
    }

    Nx = x_i.size();
    Ny = y_j.size();

    matmap_mesh = std::vector<std::vector<int>> (x_i.size(), std::vector<int>(y_j.size(), 0));

    for(int i=0; i<x_i.size(); i++)
    {
        for(int j=0; j<y_j.size(); j++)
        {
            matmap_mesh[i][j] = input_object.matmap_cell[cell_i[i]][cell_j[j]]-1;
        }
    }

}

double geometry_class::sigma_t(int i, int j, int g) const
{
    int material = matmap_mesh[i][j];
    return materials[material].sigma_t[g];
}

double geometry_class::nu_sigmaf(int i, int j, int g) const
{
    int material = matmap_mesh[i][j];
    return materials[material].nu_sigma_f[g];
}

double geometry_class::chi(int i, int j, int g) const
{
    int material = matmap_mesh[i][j];
    return materials[material].chi[g];
}

double geometry_class::sigma_s(int i, int j, int to_g, int from_g) const
{
    int material = matmap_mesh[i][j];
    return materials[material].sigma_s[to_g][from_g];
}

void geometry_class::calculate_fission_density_g(std::vector<std::vector<std::vector<double>>>  const &flux, std::vector<std::vector<std::vector<double>>> &fission_density_g, double const &keff)
{
    double total_fission_density;
        for(int i=0; i<Nx; i++)
        {
            for (int j=0; j<Ny; j++)
            {
                total_fission_density = 0.0;
                for(int g=0; g<groups; g++)
                {
                    total_fission_density = total_fission_density + flux[g][i][j] * nu_sigmaf(i,j,g);

                }

                for(int g=0; g<groups; g++)
                {
                    fission_density_g[g][i][j] = 1.0/keff*chi(i,j,g) * total_fission_density;

                }

            }

        }
}
