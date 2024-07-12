#include <iostream>
#include<vector>
#include<cmath>

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
    

    public:
    double sigma_t= 1;
    double sigma_s = 0.7;
    double nu_sigma_f = 0.39;

    int Nx = 10, Ny = 10;

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


class boundary_flux
{
    public:
    vector<vector<vector<double>>> bottom, right, top, left ;
    
};

class sweep_direction_class
{
    //this class provides the set of data required for setting direction for sweeps. there can be left to right, bottom to top and any order for both the axes
    //for example sweep('L', 'B', geometry) will sweep from left to right and bottom to top
    public:
    int x_start, x_end, y_start, y_end , x_dir, y_dir;
    int Nx, Ny, V_reflect, H_reflect;

    sweep_direction_class(char x_first, char y_first, geometry_data* geometry)
    {
        Nx = geometry->Nx;
        Ny = geometry->Ny;

        
        //x axis
        if (x_first = 'L' )
        {
            x_start = 0;
            x_end = Nx;
            x_dir = 1;

            if (y_first = 'B' )
            {
                y_start = 0;
                y_end =0;
                y_dir = 1;
                V_reflect = 1;
                H_reflect = 3;
            }
            if (y_first = 'T' )
            {
                y_start = geometry->Ny;
                y_end = 0;
                y_dir = -1;
                V_reflect = 2;
                H_reflect = 0;
            }
        }
        if (x_first = 'R' )
        {
            x_start = geometry->Nx;
            x_end = 0;
            x_dir = -1;

            if (y_first = 'B' )
            {
                y_start = 0;
                y_end =0;
                y_dir = 1;
                V_reflect = 0;
                H_reflect = 2;
            }
            if (y_first = 'T' )
            {
                y_start = geometry->Ny;
                y_end = 0;
                y_dir = -1;
                V_reflect = 3;
                H_reflect = 1;
            }
        }

        //y axis

    }
};



class sweeper_class{

    public:
    geometry_data geometry;
    angular angle;
    vector<vector<double>> psi_ij;
    vector<vector<double>> flux_ij;
    vector<vector<double>> Q;
    vector<vector<vector<double>>> psi_boun_B, psi_boun_R, psi_boun_T, psi_boun_L ;
    sweeper_class(geometry_data geometry1, angular angle1)
    {
        geometry = geometry1;
        angle = angle1;
        vector<vector<double>> psi_ij(2*geometry.Nx+1, vector<double> (2*geometry.Ny+1,0));
    }

    int transport_sweep(vector<vector<double>> q, vector<vector<double>> flux)
    {
        flux_ij = flux;
        Q = q;
        int n;
        for (n = 0; n< (angle.total_num); n++);
        {
            sweep_direction_class sweep_dir('L', 'B', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_direction_class sweep_dir('R', 'B', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_direction_class sweep_dir('L', 'T', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_direction_class sweep_dir('R', 'T', &geometry);
            sweep_part(sweep_dir, n, Q);
        }
        return 0;
    }


    int sweep_part(sweep_direction_class sweep_dir, int n, vector<vector<double>> Q)
    {
        int cell_x, cell_y;
        for(int i=sweep_dir.x_start; 0<=i && i<sweep_dir.Nx; i=i+sweep_dir.x_dir)
        {
            cell_x = 2*i+1;
            for(int j=sweep_dir.y_start; 0<=j && j<sweep_dir.Ny;j+= sweep_dir.y_dir)
            {
                cell_y = 2*j+1;
                psi_ij[cell_x][cell_y] = (geometry.sigma_t + 2*angle.mu[n]/geometry.del_x + 2*angle.mu[n]/geometry.del_y)*(2*angle.mu[n]/geometry.del_x * psi_ij[cell_x-sweep_dir.x_dir][cell_y] + 2*angle.mu[n]/geometry.del_y*psi_ij[cell_x][cell_y-sweep_dir.y_dir]+1/(4*M_PI)*Q[cell_x][cell_y]) ;
                psi_ij[cell_x+sweep_dir.x_dir][cell_y] = 2*psi_ij[cell_x][cell_y] - psi_ij[cell_x-sweep_dir.x_dir][cell_y];
                psi_ij[cell_x][cell_y+sweep_dir.y_dir] = 2*psi_ij[cell_x][cell_y] - psi_ij[cell_x][cell_y-sweep_dir.y_dir];
                flux_ij[i][j] = flux_ij[i][j] + (1/4)*angle.w[n]*psi_ij[cell_x][cell_y];
            }
            
        }
    }
};



int main()
{

    //angular discretization
    angular angle;
    geometry_data geometry;
    vector<vector<double>> flux(geometry.Nx,vector<double> (geometry.Ny,0));
    vector<vector<double>> flux_old(geometry.Nx,vector<double> (geometry.Ny,0));
    vector<vector<double>> q(geometry.Nx,vector<double> (geometry.Ny,0));

    



    

}
