#include <iostream>
#include<vector>
#include<cmath>
#include<iomanip>

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
        vector<int>sequence1{1,
                            2, 2, 
                            3, 5, 3, 
                            4, 6, 6, 4, 
                            4, 7, 8, 7, 4, 
                            3, 6, 8, 8, 6, 3, 
                            2, 5, 6, 7, 6, 5, 2, 
                            1, 2, 3, 4, 4, 3, 2, 1 };
        
        sequence = sequence1;
        for (int i=0; i<size(sequence);i++)
            sequence[i]=sequence[i]-1;
        
        total_num = size(sequence);

        vector<double>mu1{   0.1389568, 
                            0.3922893,
                            0.5370966, 
                            0.6504264, 
                            0.7467506, 
                            0.8319966, 
                            0.9092855, 
                            0.9805009};
        mu = mu1;
        vector<double>w1={   0.0489872, 
                            0.0413296, 
                            0.0212326, 
                            0.0256207, 
                            0.0360486, 
                            0.0144589, 
                            0.0344958,
                            0.0085179 };
        w=w1;
    }    
};

class geometry_data
{
    double X = 4.0, Y=4.0;
    

    public:
    double sigma_t= 1;
    double sigma_s = 0.7;
    double nu_sigma_f = 0.39;
    int albedo[4] = {0,0,0,0};

    int Nx = 8, Ny = 10;

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


class boundary_flux_class
{
    public:
    vector<vector<vector<vector<double>>>> bottom_top, left_right;
    boundary_flux_class(int xpoints,int ypoints, int angles)
    {
        vector<vector<vector<vector<double>>>> bottom_top1(2, vector<vector<vector<double>>>(4, vector<vector<double>> (xpoints, vector<double> (angles, 0))));
        bottom_top = bottom_top1;

        vector<vector<vector<vector<double>>>> left_right1(2, vector<vector<vector<double>>>(4, vector<vector<double>> (ypoints, vector<double> (angles, 0))));
        left_right = left_right1;
    }
    
};

class sweep_direction_class
{
    //this class provides the set of data required for setting direction for sweeps. there can be left to right, bottom to top and any order for both the axes
    //for example sweep('L', 'B', geometry) will sweep from left to right and bottom to top
    public:
    int x_start, x_end, y_start, y_end , x_dir, y_dir, V_boun_start, H_boun_start;
    int Nx, Ny, H_reflect, V_reflect, H_in, V_in;

    void set_sweep_direction(char x_first, char y_first, geometry_data* geometry)
    {
        Nx = geometry->Nx;
        Ny = geometry->Ny;
        int L=0, R=1;
        int B=0, T=1;
        
        //x axis
        if (x_first == 'L' )
        {
            x_start = 0;
            x_end = Nx-1;
            x_dir = 1;
            V_boun_start = L;

            if (y_first == 'B' )
            {
                y_start = 0;
                y_end =Ny-1;
                y_dir = 1;
                H_in = 0;
                H_reflect = 3;
                V_in = 0;
                V_reflect = 1;
                H_boun_start = B;

            }
            if (y_first == 'T' )
            {
                y_start = geometry->Ny-1;
                y_end = 0;
                y_dir = -1;
                H_in = 3;
                H_reflect = 0;
                V_in = 3;
                V_reflect = 2;
                H_boun_start = T;
            }
        }
        if (x_first == 'R' )
        {
            x_start = geometry->Nx-1;
            x_end = 0;
            x_dir = -1;
            V_boun_start = R;

            if (y_first == 'B' )
            {
                y_start = 0;
                y_end =Ny-1;
                y_dir = 1;
                H_in = 1;
                H_reflect = 2;
                V_in = 1;
                V_reflect = 0;
                H_boun_start = B;
            }
            if (y_first == 'T' )
            {
                y_start = geometry->Ny-1;
                y_end = 0;
                y_dir = -1;
                H_in = 2;
                H_reflect = 1;
                V_in = 2;
                V_reflect = 3;
                H_boun_start = T;
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
    vector<vector<vector<vector<double>>>> bottom_top, left_right;
    
    sweeper_class(geometry_data geometry1, angular angle1)
    {
        geometry = geometry1;
        angle = angle1;
        vector<vector<double>> psi_ij1(2*geometry.Nx+1, vector<double> (2*geometry.Ny+1,0));
        psi_ij = psi_ij1;

        
        vector<vector<vector<vector<double>>>> bottom_top1(2, vector<vector<vector<double>>>(4, vector<vector<double>> (geometry.Nx, vector<double> (angle.total_num, 0))));
        bottom_top = bottom_top1;

        vector<vector<vector<vector<double>>>> left_right1(2, vector<vector<vector<double>>>(4, vector<vector<double>> (geometry.Ny, vector<double> (angle.total_num, 0))));
        left_right = left_right1;
    }

    vector<vector<double>> transport_sweep(vector<vector<double>> q, vector<vector<double>> flux)
    {
        flux_ij = flux;
        Q = q;
        
        sweep_direction_class sweep_dir;
        
        for (int n = 0; n<angle.total_num; n++)
        {
            sweep_dir.set_sweep_direction('L', 'B', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_dir.set_sweep_direction('R', 'B', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_dir.set_sweep_direction('L', 'T', &geometry);
            sweep_part(sweep_dir, n, Q);

            sweep_dir.set_sweep_direction('R', 'T', &geometry);
            sweep_part(sweep_dir, n, Q);
        }
        return flux_ij;
    }


    void sweep_part(sweep_direction_class sweep_dir, int n, vector<vector<double>> Q)
    {
        
        int cell_x, cell_y;
        for(int i=sweep_dir.x_start; 0<=i && i<sweep_dir.Nx; i=i+sweep_dir.x_dir)
        {   
            cell_x = 2*i+1;
            for(int j=sweep_dir.y_start; 0<=j && j<sweep_dir.Ny;j+= sweep_dir.y_dir)
            {
                cell_y = 2*j+1;
                if(i==sweep_dir.x_start) //vertical boundary
                {
                    psi_ij[cell_x-sweep_dir.x_dir][cell_y] = left_right[sweep_dir.V_boun_start][sweep_dir.V_in][j][n];
                }
                if(j==sweep_dir.y_start) //horizontal boundary
                {
                    psi_ij[cell_x][cell_y - sweep_dir.y_dir] = bottom_top[sweep_dir.H_boun_start][sweep_dir.H_in][i][n];
                }
                psi_ij[cell_x][cell_y] = 1.0/(geometry.sigma_t + 2.0*angle.mu[angle.sequence[n]]/geometry.del_x + 2.0*angle.mu[angle.sequence[n]]/geometry.del_y)*(2.0*angle.mu[angle.sequence[n]]/geometry.del_x * psi_ij[cell_x-sweep_dir.x_dir][cell_y] + 2.0*angle.mu[angle.sequence[n]]/geometry.del_y*psi_ij[cell_x][cell_y-sweep_dir.y_dir]+1.0/(4.0*M_PI)*Q[i][j]) ;
                psi_ij[cell_x+sweep_dir.x_dir][cell_y] = 2.0*psi_ij[cell_x][cell_y] - psi_ij[cell_x-sweep_dir.x_dir][cell_y];
                psi_ij[cell_x][cell_y+sweep_dir.y_dir] = 2.0*psi_ij[cell_x][cell_y] - psi_ij[cell_x][cell_y-sweep_dir.y_dir];
                
                flux_ij[i][j] = flux_ij[i][j] + (1.0/4.0)*angle.w[angle.sequence[n]]*psi_ij[cell_x][cell_y];

                if(i==sweep_dir.x_end) //vertical boundary
                {
                    left_right[abs(1-sweep_dir.V_boun_start)][sweep_dir.V_reflect][j][n] = 0*psi_ij[cell_x+sweep_dir.x_dir][cell_y];
                }
                if(j==sweep_dir.y_end) //horizontal boundary
                {
                    bottom_top[abs(1-sweep_dir.H_boun_start)][sweep_dir.H_reflect][i][n] = 0*psi_ij[cell_x][cell_y+sweep_dir.y_dir];
                }                 
            


            }

        }
        
    }
};



int main()
{

    //angular discretization
    angular angle;
    geometry_data geometry;
    
    sweeper_class sweep_object(geometry, angle);
    vector<vector<double>> flux_old(geometry.Nx, vector<double> (geometry.Ny,1));
    vector<vector<double>> q(geometry.Nx, vector<double> (geometry.Ny,0));
    vector<vector<double>> fission_source = q, fission_source_new;
    for(int i=0; i<geometry.Nx; i++)
    {
        for (int j=0; j<geometry.Ny; j++)
        {
            fission_source[i][j] = flux_old[i][j] * geometry.nu_sigma_f;
            q[i][j] = (fission_source[i][j]+flux_old[i][j]*geometry.sigma_s); 
        }
        
    }

    sweep_object.transport_sweep(q, flux_old);
    for(int i=0; i<geometry.Nx; i++)
    {
        for (int j=0; j<geometry.Ny; j++)
        {
            fission_source[i][j] = sweep_object.flux_ij[i][j] * geometry.nu_sigma_f;
            q[i][j] = (fission_source[i][j]+sweep_object.flux_ij[i][j]*geometry.sigma_s); 
        }
        
    }

    sweep_object.transport_sweep(q, sweep_object.flux_ij);



    for(int i=0; i<geometry.Nx; i++)
    {
        for(int j=0; j<geometry.Ny; j++)
        {
            cout<<sweep_object.flux_ij[i][j]<<" ";
        }
        cout<<endl;
    }





}
