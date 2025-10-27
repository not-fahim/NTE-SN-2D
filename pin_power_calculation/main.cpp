#include <iostream>
#include<vector>
#include <cctype>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<sstream>
#include "geometry_class.hpp"
#include "geo_mat_parser.hpp"
#include "pin_power.hpp"

using namespace std;
int main ()

{
    //read input file
    input_class input_object = read_input_file();
    //create geometry object
    geometry_class geometry(input_object);
    cout<< geometry.X << "  "<< geometry.Y <<endl;
}

