#ifndef GEO_MAT_PARSER
#define GEO_MAT_PARSER
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <iomanip>

struct cell_class
{   
    double cell_dim;
    int mesh_num;
};

struct mat_class
{
    std::string name;
    std::vector<double> sigma_t, nu_sigma_f, chi;
    std::vector<std::vector<double>> sigma_s;
};

struct input_class
{
    std::string name, xsfile;
    std::vector<cell_class> cellX;
    std::vector<cell_class> cellY;
    std::vector< std::vector<int>> matmap_cell;
    std::vector<double> BC;
    int S_n;
    std::vector<mat_class> mat_vector;
};

input_class read_input_file();
void parseCoordinates( std::stringstream& ss,  std::vector<cell_class>& targetVector);
std::vector<mat_class> read_material_input(std::string filename);
void print_input_data(const input_class& data);

#endif

