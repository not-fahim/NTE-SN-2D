#include "pin_power.hpp"
#include <iostream> // For std::cerr (logging errors)
#include <fstream>
#include <sstream>
#include <iomanip> // For std::setprecision

pin_class::pin_class(geometry_class const &geometry, input_class const &input_object)
{
    int pin_num_x = input_object.pin_dim[0];
    int pin_num_y = input_object.pin_dim[1];
    pin_pitch = input_object.pin_pitch;

    pin_volume_inp = pin_pitch * pin_pitch;

    pin_xi_plush_half = range(pin_pitch, geometry.X, pin_num_x);
    pin_xi_minush_half = range(0, geometry.X-pin_pitch, pin_num_x);
    pin_xi = cell_midpoints(pin_xi_minush_half, pin_xi_plush_half);

    pin_yj_plush_half = range(pin_pitch, geometry.Y, pin_num_y);
    pin_yj_minush_half = range(0, geometry.Y-pin_pitch, pin_num_y);
    pin_yj = cell_midpoints(pin_yj_minush_half, pin_yj_plush_half);

    pin_fission_rate = std::vector<std::vector<std::vector<double>>> (geometry.groups, std::vector<std::vector<double>>(pin_num_x, std::vector<double>(pin_num_y, 0.0)));
    pin_flux = std::vector<std::vector<std::vector<double>>> (geometry.groups, std::vector<std::vector<double>>(pin_num_x, std::vector<double>(pin_num_y, 0.0)));
    pin_volume_mesh = std::vector<std::vector<double>> (pin_num_x, std::vector<double>(pin_num_y, 0.0));

}

