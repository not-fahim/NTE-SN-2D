#ifndef PIN_POWER_CALCULATION_HPP
#define PIN_POWER_CALCULATION_HPP
#include <vector>
#include <string>

struct doutb_flux_output_class
{
    std::string name;
    int groups;
    bool pin_calc = false;
    std::vector<int> pin_num;
    double pin_pitch=0.0;
    std::vector<std::vector<std::vector<double>>> mesh_flux; // [group][mesh_x][mesh_y]
    std::vector<double> mesh_xi, mesh_yj;
    std::vector<double> mesh_delxi, mesh_delyj;
    std::vector<std::vector<double>> mesh_fission_density;
};

doutb_flux_output_class read_dotb_flux_ouput(std::string filename);

class pin_data_class
{
    public:
    
    // Pin geometry info,  will come from theinput file
    double pin_pitch;
    double pin_volume;
    int pin_num_x, pin_num_y;
    int groups;

    //pin coordinates and boundaries will be calculated from pin pitch and number of pins
    std::vector<double> pin_xi, pin_yj; // pin midpoints in x and y
    std::vector<double> pin_xi_minus_half, pin_yj_minus_half; // pin boundaries in x and y
    std::vector<double> pin_xi_plus_half, pin_yj_plus_half;// pin boundaries in x and y

    //pin data storage
    std::vector<std::vector<std::vector<double>>>  pin_flux; // [group][pin_x][pin_y]
    std::vector<std::vector<double>> pin_volume_mesh, pin_fission_rate; //fission rate [pin_x][pin_y]

    // mesh info from the solver, comes from the input file
    std::vector<double> mesh_xi, mesh_yj; // solver mesh cell midpoints in x and y
    std::vector<double> mesh_delxi, mesh_delyj; // solver mesh cell widths in x and y
    double total_fission_rate;
    std::vector<std::vector<std::vector<double>>> mesh_flux; // [group][mesh_x][mesh_y]
    std::vector<std::vector<double>> mesh_fission_density;

    pin_data_class(doutb_flux_output_class const &dotb);
    void calculate_pin_avg_group_flux();
    void calculate_pin_fission_rate_from_mesh();



    
};

#endif