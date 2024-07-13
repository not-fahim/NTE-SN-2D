#include <iostream>
#include<vector>

using namespace std;
vector<double> range(double min, double max, size_t N) 
{
    vector<double> range;
    double delta = (max-min)/double(N-1);
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    
    return range;
}

vector<double> cell_midpoints(vector<double> x_i_minus_half, vector<double> x_i_plus_half)
{
    vector<double> cell_midpoints;
    for(int i=0; i<size(x_i_minus_half); i++)
    {
        cell_midpoints.push_back((x_i_minus_half[i] + x_i_plus_half[i])/2);
    }
    return cell_midpoints;
}

class geometry_data
{
    double X = 4.0, Y=4.0;
    

    public:
    double sigma_t= 1;
    double sigma_s = 0.7;
    double nu_sigma_f = 0.39;
    int albedo[4] = {0,0,0,0};

    int Nx = 9, Ny = 10;

    double del_x = X/ double(Nx), del_y = Y/double(Ny);
    vector<double> x_i_minus_half, x_i_plus_half, x_i ;
    vector<double> y_j_minus_half, y_j_plus_half, y_j;

    

    geometry_data()
    {
        x_i_minus_half = range(0,X-del_x,Nx);
        x_i_plus_half = range(del_x, X, Nx);
        x_i = cell_midpoints(x_i_minus_half, x_i_plus_half);

        y_j_minus_half = range(0,Y-del_y,Ny);
        y_j_plus_half = range(del_y, Y, Ny);
        y_j = cell_midpoints(y_j_minus_half,x_i_plus_half);   



    }

};

int main()
{
    geometry_data geometry;
    for(int i=0; i<size(geometry.x_i); i++) {
        
        cout<<geometry.x_i[i]<<'\t';
    }
    

}