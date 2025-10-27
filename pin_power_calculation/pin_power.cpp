#include "pin_power.hpp"
#include <iostream> // For std::cerr (logging errors)
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip> // For std::setprecision

std::vector<double> range(double min, double max, size_t N)
{
    //returns a vector of N equally spaced doubles including min  max
    std::vector<double> range;
    double delta = (max-min)/double((N>1)? N-1 : 1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    return range;
}

std::vector<double> cell_midpoints(std::vector<double> x_i_minus_half, std::vector<double> x_i_plus_half)
{
    std::vector<double> cell_midpoints;
    for(int i=0; i<size(x_i_minus_half); i++)
    {
        cell_midpoints.push_back((x_i_minus_half[i] + x_i_plus_half[i])/2);
    }
    return cell_midpoints;
}


doutb_flux_output_class read_dotb_flux_ouput(std::string filename)
{
    bool pincalc = false;
    doutb_flux_output_class dotb_output;
    int groups;

    std::vector<std::vector<double>> group_flux_temp;

    std::ifstream inputFile(filename); 
    
    if (!inputFile.is_open()) 
    {
         std::cerr << "Error: Could not open:"<<filename <<  std::endl;
        dotb_output.name = "nan"; //works as flag for the error
        return dotb_output;
    }

    std::string line;
    // one line at a time
    while ( getline(inputFile, line)) 
    {
        // Skip any empty lines
        if (line.empty()) 
        {
            continue; 
        }
        std::stringstream linestream(line);
        std::string keyword;
        linestream >> keyword; //first word to identify the line's purpose

        if(keyword[0] == '%') //comment line
        {
            continue;
        }
        // the line based on its keyword
        else if(keyword=="name")
        {
            std::string name;
            if(linestream >> name)
                dotb_output.name = name;
        }

        else if(keyword=="pincalc")
        {
            std::string value;
            if(linestream >> value)
            {
                if(value == "true")
                {
                    pincalc = true;
                    dotb_output.pin_calc = true;
                }
            }
        }

        else if(keyword == "pinpitch" && pincalc==true) 
        {
            double value;
            if(linestream >> value)
                dotb_output.pin_pitch = value;
        }

        else if(keyword == "pins" && pincalc==true) 
        {
            int value;
            while(linestream >> value)
                dotb_output.pin_num.push_back(value);
        }

        else if (keyword == "cellx") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_xi.push_back(value);
        } 
        else if (keyword == "celly") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_yj.push_back(value);
        } 
        else if (keyword == "delxi") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_delxi.push_back(value);
        }
        else if (keyword == "delyj") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_delyj.push_back(value);
        }
        else if (keyword == "groups") 
        {
            int value;
            if(linestream >> value)
                dotb_output.groups = value;
        }
        else if(keyword == "fiss_den")
        {
            std::vector<std::vector<double>> fission_density_temp;   
            std::string matrixLine;
            while ( getline(inputFile, matrixLine) && !matrixLine.empty()) 
            {
                std::stringstream matrixStream(matrixLine);
                std::vector<double> row;
                double value;
                while (matrixStream >> value) 
                {
                    row.push_back(value);
                }
                if (!row.empty()) 
                {
                    fission_density_temp.push_back(row);
                }
                
            }
            dotb_output.mesh_fission_density = fission_density_temp;
        }
        for(int g=0; g<groups; g++)
        {
            if (keyword == "phi"+std::to_string(g+1)) 
            {
                group_flux_temp.clear();
                std::string matrixLine;
                while ( getline(inputFile, matrixLine) && !matrixLine.empty()) 
                {
                    std::stringstream matrixStream(matrixLine);
                    std::vector<double> row;
                    double value;
                    while (matrixStream >> value) 
                    {
                        row.push_back(value);
                    }
                    if (!row.empty()) 
                    {
                        group_flux_temp.push_back(row);
                    }
                }
                dotb_output.mesh_flux.push_back(group_flux_temp);
            }
        }
        
    }
    inputFile.close();
    return dotb_output;

}

pin_data_class::pin_data_class(doutb_flux_output_class const &dotb)
{
    groups = dotb.groups;
    pin_num_x = dotb.pin_num[0];
    pin_num_y = dotb.pin_num[1];
    pin_pitch = dotb.pin_pitch;
    pin_volume = pin_pitch * pin_pitch;

    mesh_xi = dotb.mesh_xi;
    mesh_yj = dotb.mesh_yj;
    mesh_delxi = dotb.mesh_delxi;
    mesh_delyj = dotb.mesh_delyj;
    mesh_flux = dotb.mesh_flux;
    total_fission_rate = 0.0;
    mesh_fission_density = dotb.mesh_fission_density;

    double X = pin_pitch*double(pin_num_x);
    double Y = pin_pitch*double(pin_num_y);

    pin_xi_plus_half = range(pin_pitch, X, pin_num_x);
    pin_xi_minus_half = range(0, X-pin_pitch, pin_num_x);
    pin_xi = cell_midpoints(pin_xi_minus_half, pin_xi_plus_half);

    pin_yj_plus_half = range(pin_pitch, Y, pin_num_y);
    pin_yj_minus_half = range(0, Y-pin_pitch, pin_num_y);
    pin_yj = cell_midpoints(pin_yj_minus_half, pin_yj_plus_half);

    pin_fission_rate = std::vector<std::vector<double>>(pin_num_x, std::vector<double>(pin_num_y, 0.0));

    pin_flux = std::vector<std::vector<std::vector<double>>> (groups, std::vector<std::vector<double>>(pin_num_x, std::vector<double>(pin_num_y, 0.0)));
    pin_volume_mesh = std::vector<std::vector<double>> (pin_num_x, std::vector<double>(pin_num_y, 0.0));

}

doutb_flux_output_class read_dotb_flux_ouput(std::string filename)
{
    bool pincalc = false;
    doutb_flux_output_class dotb_output;
    int groups;

    std::vector<std::vector<double>> group_flux_temp;

    std::ifstream inputFile(filename); 
    
    if (!inputFile.is_open()) 
    {
         std::cerr << "Error: Could not open:"<<filename <<  std::endl;
        dotb_output.name = "nan"; //works as flag for the error
        return dotb_output;
    }

    std::string line;
    // one line at a time
    while ( getline(inputFile, line)) 
    {
        // Skip any empty lines
        if (line.empty()) 
        {
            continue; 
        }
        std::stringstream linestream(line);
        std::string keyword;
        linestream >> keyword; //first word to identify the line's purpose

        if(keyword[0] == '%') //comment line
        {
            continue;
        }
        // the line based on its keyword
        else if(keyword=="name")
        {
            std::string name;
            if(linestream >> name)
                dotb_output.name = name;
        }

        else if(keyword=="pincalc")
        {
            std::string value;
            if(linestream >> value)
            {
                if(value == "true")
                {
                    pincalc = true;
                    dotb_output.pin_calc = true;
                }
            }
        }

        else if(keyword == "pinpitch" && pincalc==true) 
        {
            double value;
            if(linestream >> value)
                dotb_output.pin_pitch = value;
        }

        else if(keyword == "pins" && pincalc==true) 
        {
            int value;
            while(linestream >> value)
                dotb_output.pin_num.push_back(value);
        }

        else if (keyword == "cellx") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_xi.push_back(value);
        } 
        else if (keyword == "celly") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_yj.push_back(value);
        } 
        else if (keyword == "delxi") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_delxi.push_back(value);
        }
        else if (keyword == "delyj") 
        {
            double value;
            while(linestream >> value)
                dotb_output.mesh_delyj.push_back(value);
        }
        else if (keyword == "groups") 
        {
            int value;
            if(linestream >> value)
                dotb_output.groups = value;
        }
        for(int g=0; g<groups; g++)
        {
            if (keyword == "phi"+std::to_string(g+1)) 
            {
                group_flux_temp.clear();
                std::string matrixLine;
                while ( getline(inputFile, matrixLine) && !matrixLine.empty()) 
                {
                    std::stringstream matrixStream(matrixLine);
                    std::vector<double> row;
                    double value;
                    while (matrixStream >> value) 
                    {
                        row.push_back(value);
                    }
                    if (!row.empty()) 
                    {
                        group_flux_temp.push_back(row);
                    }
                }
                dotb_output.mesh_flux.push_back(group_flux_temp);
            }
        }
        
    }
    inputFile.close();
    return dotb_output;

}

