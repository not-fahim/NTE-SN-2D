#include <iostream>
#include <string>
#include <iomanip> // For std::setprecision
#include <fstream>
#include "pin_power.hpp"
#include <sstream>

int main(char argc, char* argv[])
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <dotb_flux_output_filename>" << std::endl;
        return 1;
    }

    std::string dotb_filename = argv[1];
    pin_data_class pin_data(dotb_filename);
    std::cout << "Pin power calculation started..." << std::endl;
    pin_data.calculate_pin_power();
    std::cout << "Pin power calculation done, writing output file..." << std::endl;
    pin_data.write_pin_power_output();

    return 0;
}