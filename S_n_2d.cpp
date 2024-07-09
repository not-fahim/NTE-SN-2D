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

class angular
{
    public:
    vector<int> sequence;
    vector<double> mu;
    vector<double> w;
    int total_num;
    angular()
    {
        vector<int>sequence{1,
                            2, 2, 
                            3, 5, 3, 
                            4, 6, 6, 4, 
                            4, 7, 8, 7, 4, 
                            3, 6, 8, 8, 6, 3, 
                            2, 5, 6, 7, 6, 5, 2, 
                            1, 2, 3, 4, 4, 3, 2, 1 };
        for (int i=0; i<size(sequence);i++)
            sequence[i]=sequence[i]-1;
        total_num = size(sequence);
        vector<double>mu{   0.1389568, 
                            0.3922893,
                            0.5370966, 
                            0.6504264, 
                            0.7467506, 
                            0.8319966, 
                            0.9092855, 
                            0.9805009};

        vector<double>w={   0.0489872, 
                            0.0413296, 
                            0.0212326, 
                            0.0256207, 
                            0.0360486, 
                            0.0144589, 
                            0.0344958,
                            0.0085179 };
    }    
};

class geometry_data
{
    double X = 4.0, Y=4.0;
    int Nx = 10, Ny = 10;

    public:

    double del_x = X/ double(Nx), del_y = Y/double(Ny);
    vector<double> x_i_minus_half, x_i_plus_half, x_i ;
    vector<double> y_j_minus_half, y_j_plus_half, y_j;

    geometry_data()
    {
        vector<double> x_i_minus_half = range(0,X-del_x,Nx);
        vector<double> x_i_plus_half = range(del_x, X, Nx);
        vector<double> x_i = cell_midpoints(x_i_minus_half, x_i_plus_half);

        vector<double> y_j_minus_half = range(0,Y-del_y,Ny);
        vector<double> y_j_plus_half = range(del_y, Y, Ny);
        vector<double> y_j = cell_midpoints(y_j_minus_half,x_i_plus_half);   
    }

};

class flux_and_boundary_flux_data
{
    public:
    double* flux, boundary_flux;

    flux_and_boundary_flux_data(int Nx, int Ny)
    {
        double flux[Nx][Ny];
        double boundary_flux[Nx][Ny];
    }
};

int main()
{
    double sigma_t= 1;
    double sigma_s = 0.7;
    double nu_sigma_f = 0.39;

    //angular discretization
    angular angle;
    geometry_data geometry;


    

}
