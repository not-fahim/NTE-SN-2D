#ifndef GEO_MAT_PARSER
#define GEO_MAT_PARSER
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>

struct cell_class
{   
    double cell_dim;
    int mesh_num;
};

struct mat_class
{
    std::string name;
    std::vector<double> tot, nuf, chi, sca;
};

struct input_class
{
    std::string name;
    std::vector<cell_class> cellX;
    std::vector<cell_class> cellY;
    std::vector< std::vector<int>> matmap_cell;
    std::vector<double> BC;
    int S_n;
    std::vector<mat_class> materials;
};

input_class read_input_file();
void parseCoordinates( std::stringstream& ss,  std::vector<cell_class>& targetVector);

#endif