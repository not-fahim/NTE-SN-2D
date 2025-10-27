#ifndef PIN_POWER_CALCULATION_HPP
#define PIN_POWER_CALCULATION_HPP
#include <vector>
#include <string>
#include "geometry_class.hpp"
#include "geo_mat_parser.hpp"

struct mesh_flux_data
{
    std::vector<std::vector<std::vector<double>>> mesh_flux; // [group][mesh_x][mesh_y]
    std::vector<double> mesh_xi, mesh_yj;
    std::vector<double> mesh_delxi, mesh_delyj;
};

mesh_flux_data read_mesh_flux_file(std::string filename);

class pin_data_class
{
    public:
    double pin_pitch;
    double pin_volume_inp;

    std::vector<double> pin_xi, pin_yj; // pin midpoints in x and y
    std::vector<double> pin_xi_minush_half, pin_yj_minush_half; // pin boundaries in x and y
    std::vector<double> pin_xi_plush_half, pin_yj_plush_half;// pin boundaries in x and y

    std::vector<std::vector<std::vector<double>>>  pin_flux; // [group][pin_x][pin_y]
    std::vector<std::vector<double>> pin_volume_mesh, pin_fission_rate; //fission rate [pin_x][pin_y]
    double total_fission_rate;
    void read_flux_file();
    void calculate_pin_powers(geometry_class const &geometry, input_class const &input_object);
    pin_data_class(geometry_class const &geometry, input_class const &input_object);
};

#endif