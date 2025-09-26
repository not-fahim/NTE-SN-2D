#ifndef PLOTTING_HPP
#define PLOTTING_HPP
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "geometry_class.hpp"

void write_matlab_flux_script(const std::vector<std::vector<std::vector<double>>>& flux, const geometry_class& geometry, const std::string& problem_name);
void write_python_notebook_script(const std::vector<std::vector<std::vector<double>>>& flux, const geometry_class& geometry, const std::string& problem_name);

#endif