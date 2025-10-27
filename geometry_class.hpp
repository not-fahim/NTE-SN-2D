#ifndef GEOMETRY_H
#define GEOMETRY_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>
#include "geo_mat_parser.hpp"

class geometry_class
{
    public:

    double X;
    double Y;
    int Nx;
    int Ny;
    int groups;

    std::vector<double> x_i_minus_half, x_i_plus_half, x_i, delx_i ;
    std::vector<double> y_j_minus_half, y_j_plus_half, y_j, dely_j;
    std::vector<std::vector<int>> matmap_mesh;
    std::vector<mat_class> materials;

    //boundary condition:
    double alpha_B, alpha_T, alpha_L, alpha_R;

    //methods:
    geometry_class(input_class const &input_object);
    double sigma_t(int i, int j, int g) const;
    double sigma_s(int i, int j, int to_g, int from_g) const;
    double nu_sigmaf(int i, int j, int g) const;
    double chi(int i,int j,int g) const;
    double sigma_f(int i,int j,int g) const;
    void  calculate_fission_density_g(std::vector<std::vector<std::vector<double>>>  const &flux, std::vector<std::vector<std::vector<double>>> &fission_density_g, double const &keff);
};

int vector3Dcopy(std::vector<std::vector<std::vector<double>>> &to_vector1, std::vector<std::vector<std::vector<double>>> const &from_vector2);
std::vector<double> range(double min, double max, size_t N);
std::vector<double> cell_midpoints(std::vector<double> x_i_minus_half, std::vector<double> x_i_plus_half);
#endif