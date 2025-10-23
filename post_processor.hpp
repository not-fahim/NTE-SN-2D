#ifndef POST_PROCESSOR_H
#define POST_PROCESSOR_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "geo_mat_parser.hpp"
#include "geometry_class.hpp"
class post_processor_class
{
    public:
    geometry_class geometry;
    std::vector<std::vector<std::vector<double>>> flux_g; //3D vector: group,x,y
    std::vector<std::vector<std::vector<double>>> fission_density_g; //3D vector: group,x,y
    double keff;

    //methods
    post_processor_class(geometry_class const &geometry_object,
                         std::vector<std::vector<std::vector<double>>> const &flux_input, 
                         input_class const &input_object);
    void output_flux_to_file(std::string filename_prefix);
    void output_fission_density_to_file(std::string filename_prefix);
};
#endif