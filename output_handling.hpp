#ifndef OUTPUT_WRITER_HPP
#define OUTPUT_WRITER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "geometry_class.hpp" // Contains geometry_class definition
#include "geo_mat_parser.hpp" // Contains input_class definition

/**
 * @brief Writes the final scalar flux and geometry information to an output file.
 *
 * The output file will be named "flux" + input.name + ".out".
 * The file contains:
 * - A comment header describing data organization.
 * - Cell center coordinates in x (cellx) and y (celly).
 * - The refinement level.
 * - The scalar flux for each group (phi1, phi2, ...), printed as a 2D matrix.
 * In the flux matrix, columns represent the x-index (i) and rows represent the y-index (j).
 *
 * @param input     The input_class object, containing the problem name and refinement level.
 * @param geometry  The geometry_class object, containing mesh info (Nx, Ny, groups, x_i, y_j).
 * @param flux      The 3D vector of scalar flux, assumed to be indexed as [g][i][j].
 */
void write_flux_output(
    const input_class& input,
    const geometry_class& geometry,
    const std::vector<std::vector<std::vector<double>>>& flux
);
void print_banner(std::ofstream &outputfile);
#endif // OUTPUT_WRITER_HPP
