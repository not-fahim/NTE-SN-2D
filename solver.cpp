#include <iostream>
#include<vector>
#include<cmath>
#include<iomanip>
#include<fstream>
#include "UAL-2-0-SL.cpp"


double intol = 10e-10;
double outtol=10e-10;
double max_it = 5000;
using namespace std;

template <typename t>
int set_zero(vector<vector<t>> &vector1)
{
    for (int i=0; i<size(vector1); i++)
    {
        for(int j=0; j<size(vector1[i]); j++)
        {
            vector1[i][j] = 0.0;
        }
    }

    return 0;
}

template <typename t>
int set_zero(vector<vector<vector<t>>> &vector1)
{
    for (int i=0; i<size(vector1); i++)
    {
        for(int j=0; j<size(vector1[i]); j++)
        {
            for(int k=0; k<size(vector1[i][j]); k++)
            {
                vector1[i][j][k] = 0.0;
            }
        }
    }

    return 0;
}

template <typename t>
t sum(vector<vector<t>> const &vector1)
{
    t sum=0.0;
    for (int i=0; i<size(vector1); i++)
    {
        for(int j=0; j<size(vector1[0]); j++)
        {
            sum = sum+vector1[i][j];
        }
    }

    return sum;
}


double abs_max_fraction_error(vector<vector<vector<double>>> const &vector1, vector<vector<vector<double>>> const &vector2)
{

        double avg = 0.;
        double abs_max=0;
        double a;

        avg = avg / (size(vector2) * size(vector2[0]));

        for (int i=0; i< size(vector1); i++)
        {
            for (int j=0; j< size(vector1[i]); j++)
            {
                for(int k=0; k<size(vector1[i][j]); k++)
                {
                    if(i==0 && j==0 && k==0) abs_max = (vector1[i][j][k] - vector2[i][j][k]);
                
                    a=abs(vector1[i][j][k] - vector2[i][j][k]);

                    if(a>abs_max) abs_max = a;
                } 
            }
        }
        return abs_max/sizeof(vector1) * sizeof(double);

}

template <typename t>
bool is_flux_converged(vector<vector<vector<t>>> const &vector1, vector<vector<vector<t>>> const &vector2 )
{
    if(abs_max_fraction_error(vector1, vector2) < intol )
    {
        return true;
    }
    else
    {
        return false;
    }
}


double abs_max_fraction_error(vector<vector<double>> const &vector1, vector<vector<double>> const &vector2)
{

        double avg = 0.;
        double abs_max=0;
        double a;


        for (int i=0; i< size(vector1); i++)
        {
            for (int j=0; j< size(vector1[0]); j++)
            {
                if(i==0 && j==0) abs_max = (vector1[i][j] - vector2[i][j])/avg;
                
                a=abs(vector1[i][j] - vector2[i][j])/avg;

                if(a>abs_max) abs_max = a; 
            }
        }
        return abs_max;

}

int vector3Dcopy(vector<vector<vector<double>>> &to_vector1, vector<vector<vector<double>>> const &from_vector2)
{
    for(int i =0; i<to_vector1.size(); i++)
        for(int j=0; j<to_vector1[0].size(); j++)
            for(int k=0; k<to_vector1[0][0].size(); k++)
                to_vector1[i][j][k] = from_vector2[i][j][k];
    
    return 0;

}

template <typename t>
t abs_max(vector<vector<t>> const &vector1)
{
    t max = abs(vector1[0][0]);
    for(int i=0; i<size(vector1); i++)
    {
        for(int j=0; j<size(vector1[0]); j++)
        {
            if(abs(vector1[i][j]) > max )
                max = abs(vector1[i][j]);
        }

    }
    return max;
}

template <typename t>
bool is_flux_converged(vector<vector<t>> const &vector1, vector<vector<t>> const &vector2 )
{
    if(abs_max_fraction_error(vector1, vector2) < intol )
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename t>
bool is_keff_converged(t var1, t var2 )
{
    if(abs(var2-var1)/var2 < outtol )
    {
        return true;
    }
    else
    {
        return false;
    }
}

vector<double> range(double min, double max, size_t N)
{
    //returns a vector of N equally spaced doubles including min  max
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
    vector<double> mu, mu_sequence, eta_sequence;
    vector<double> w;
    int total_num;
    angular()
    {

        vector<int>
                   sequence1{ 1,
                              2, 2,
                              3, 5, 3,
                              4, 6, 6, 4,
                              4, 7, 8, 7, 4,
                              3, 6, 8, 8, 6, 3,
                              2, 5, 6, 7, 6, 5, 2,
                              1, 2, 3, 4, 4, 3, 2, 1 };


        for(int i=0; i<8 ; i++)
        {
            for(int j=0;j<=i;j++)
            {
                mu_sequence.push_back(i-j);
                eta_sequence.push_back(j);
            }
        }

        sequence = sequence1;
        for (int i=0; i<size(sequence);i++)
            sequence[i]=sequence[i]-1;

        total_num = size(sequence);

        vector<double>mu1{  0.13344572,
                            0.39119433,
                            0.53689687,
                            0.65075610,
                            0.74746822,
                            0.83302700,
                            0.91058181,
                            0.98203079};
        mu = mu1;
        vector<double>w1={  0.05415425,
                            0.03679653,
                            0.02777273,
                            0.02580284,
                            0.02494275,
                            0.01962325,
                            0.01879762,
                            0.01544801 };
        w=w1;
    }
};


class geometry_class
{
    public:

    double X;
    double Y;
    int Nx;
    int Ny;
    int groups;

    vector<double> x_i_minus_half, x_i_plus_half, x_i, delx_i ;
    vector<double> y_j_minus_half, y_j_plus_half, y_j, dely_j;


    //boundary condition:
    double alpha_B, alpha_T, alpha_L, alpha_R;

    

    geometry_class(int mesh_num, double r_c)
    {
        X= r_c;
        Y= 2;
        Nx= mesh_num;
        Ny=2;
        alpha_B = BC[0]; // bottom
        alpha_T = BC[1]; //top

        alpha_L = BC[2] ;//left
        alpha_R = BC[3];//right

        groups = sizeof(sigma_t_fuel)/sizeof(double);
        double del_x = X/ double(Nx), del_y = Y/double(Ny);

        x_i_minus_half = range(0,X-del_x,Nx);
        x_i_plus_half = range(del_x, X, Nx);
        x_i = cell_midpoints(x_i_minus_half, x_i_plus_half);


        y_j_minus_half = range(0,Y-del_y,Ny);
        y_j_plus_half = range(del_y, Y, Ny);
        y_j = cell_midpoints(y_j_minus_half,y_j_plus_half);

        cout<< "cell midpoints generated"<<endl;



        cout<<"generating cell widths"<<endl;

        for (int i = 0; i < x_i.size(); i++)
        {
            delx_i.push_back(x_i_plus_half[i] - x_i_minus_half[i]);
        }

        for (int j = 0; j < y_j.size(); j++)
        {
            dely_j.push_back(y_j_plus_half[j]- y_j_minus_half[j] );
        }

    }
    int material(int i, int j) const;

    double sigma_t(int i, int j, int g) const
    {
        return sigma_t_fuel[g];
    }

    double sigma_s(int i, int j, int to_g, int from_g) const
    {
        return sigma_s_fuel[to_g][from_g];
    }

    double nu_sigmaf(int i, int j, int g) const
    {
        return nu_sigmaf_fuel[g];

    }

    double chi(int i,int j,int g) const
    {
        return chi_fuel[g];
    }


        int  calculate_fission_density_g(vector<vector<vector<double>>>  const &flux, vector<vector<vector<double>>> &fission_density_g, double const &keff)
    {

        double total_fission_density;
        for(int i=0; i<Nx; i++)
        {
            for (int j=0; j<Ny; j++)
            {
                total_fission_density = 0.0;
                for(int g=0; g<groups; g++)
                {
                    total_fission_density = total_fission_density + flux[g][i][j] * nu_sigmaf(i,j,g);

                }

                for(int g=0; g<groups; g++)
                {
                    fission_density_g[g][i][j] = 1.0/keff*chi(i,j,g) * total_fission_density;

                }

            }

        }

        return 1;
    }

};

int geometry_class:: material(int i, int j) const
{
    return 0;
};

class Transport_sweep_class
{
    public:
    vector<vector<vector<double>>> L_boun_psi, R_boun_psi, B_boun_psi, T_boun_psi; //interim storage
    vector<vector<vector<vector<double>>>> Storage_L_boun_psi, Storage_R_boun_psi, Storage_B_boun_psi, Storage_T_boun_psi;
    double mu, eta, w;

    Transport_sweep_class(angular const &angle, geometry_class const &geometry)
    {
        //all the boundaries has two orientation: positve or negative. I am visualizing boundaries as axes and their sign is determined by the sign of tan in that quadrant
        Storage_L_boun_psi = vector<vector<vector<vector<double>>>> (geometry.groups, vector<vector<vector<double>>>(geometry.Ny, vector<vector<double>> (angle.total_num, vector<double> (2, 0.))));
        Storage_R_boun_psi = vector<vector<vector<vector<double>>>> (geometry.groups, vector<vector<vector<double>>>(geometry.Ny, vector<vector<double>> (angle.total_num, vector<double> (2, 0.))));

        Storage_B_boun_psi = vector<vector<vector<vector<double>>>> (geometry.groups, vector<vector<vector<double>>>(geometry.Nx, vector<vector<double>> (angle.total_num, vector<double> (2, 0.))));
        Storage_T_boun_psi = vector<vector<vector<vector<double>>>> (geometry.groups, vector<vector<vector<double>>>(geometry.Nx, vector<vector<double>> (angle.total_num, vector<double> (2, 0.))));
        //phi = vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny, 0.0));
        
    };


    int do_sweep (geometry_class const &geometry, angular const &angle, vector<vector<double>> const &q, int g, vector<vector<double>> &phi)
    {
        set_zero(phi);
        vector<double> psi_left_j(geometry.Ny, 0.0), psi_right_j(geometry.Ny, 0.0), psi_up_i(geometry.Nx, 0.0), psi_down_i(geometry.Nx, 0.0);
        double psi_mid;


        for(int n=0; n<angle.total_num; n++)//mu>0, eta>0, + +
        {
            mu = angle.mu[angle.mu_sequence[n]];
            eta = angle.mu[angle.eta_sequence[n]];
            w = angle.w[angle.sequence[n]];

            for(int i=0; i<geometry.Nx; i++)
            {
                for(int j=0; j<geometry.Ny; j++)
                {
                    if(i==0)
                    {
                        psi_left_j[j] = Storage_L_boun_psi[g][j][n][0];
                    }
                    else psi_left_j[j] = psi_right_j[j];

                    if(j==0)
                    {
                        psi_down_i[i] = Storage_B_boun_psi[g][i][n][0];
                    }
                    else psi_down_i[i] = psi_up_i[i];

                    psi_mid =(  2.0*mu/geometry.delx_i[i] * psi_left_j[j]
                            + 2.0*eta/geometry.dely_j[j]*psi_down_i[i]
                            + q[i][j])
                            /(
                                    geometry.sigma_t(i,j, g)
                                    + 2.0 * mu /geometry.delx_i[i]
                                    + 2.0 * eta/geometry.dely_j[j]
                                );

                    psi_right_j[j] =(2.0 * psi_mid - psi_left_j[j]);
                    
                    psi_up_i[i] = (2.0 * psi_mid - psi_down_i[i]);
                    

                    phi[i][j] = phi[i][j]+0.25*w*psi_mid;

                    if(i==(geometry.Nx-1))
                    {
                        Storage_R_boun_psi[g][j][n][1] = geometry.alpha_R*psi_right_j[j];
                    }
                    if(j==(geometry.Ny-1))
                    {
                        Storage_T_boun_psi[g][i][n][1]= geometry.alpha_T*psi_up_i[i];
                    }
                }
            }
        }
        for(int n=0; n<angle.total_num; n++)//mu<0 eta>0 - +
        {
            mu = angle.mu[angle.mu_sequence[n]];
            eta = angle.mu[angle.eta_sequence[n]];
            w = angle.w[angle.sequence[n]];

            for(int i=geometry.Nx-1; i>=0; i--)
            {
                for(int j=0; j<geometry.Ny; j++)
                {
                    if(i==(geometry.Nx-1))
                    {
                        psi_right_j[j] = Storage_R_boun_psi[g][j][n][1];
                    }
                    else psi_right_j[j] = psi_left_j[j];

                    if(j==0)
                    {
                        psi_down_i[i] = Storage_B_boun_psi[g][i][n][1];
                    }
                    else psi_down_i[i] = psi_up_i[i];

                    psi_mid = 1.0/(geometry.sigma_t(i,j,g)+ 2.0 * mu /geometry.delx_i[i] + 2.0 * eta/geometry.dely_j[j])
                              *(+2.0*mu/geometry.delx_i[i] * psi_right_j[j] + 2.0*eta/geometry.dely_j[j]*psi_down_i[i]+ q[i][j]) ;

                                
                    psi_left_j[j] =  (2.0 * psi_mid - psi_right_j[j]) ;
                    psi_up_i[i] = (2.0 * psi_mid - psi_down_i[i]);

                    phi[i][j] = phi[i][j]+0.25 * angle.w[angle.sequence[n]]*psi_mid;

                    if(i==0)
                    {
                        Storage_L_boun_psi[g][j][n][0] = geometry.alpha_L * psi_left_j[j];
                    }
                    if(j==(geometry.Ny-1))
                    {
                        Storage_T_boun_psi[g][i][n][0]= geometry.alpha_T * psi_up_i[i];
                    }
                }
            }
        }


        for(int n=0; n<angle.total_num; n++) //mu<0 eta<0 - -
        {
            mu = angle.mu[angle.mu_sequence[n]];
            eta = angle.mu[angle.eta_sequence[n]];
            w = angle.w[angle.sequence[n]];

            for(int i=geometry.Nx-1; i>=0 ; i--)
            {
                for(int j=geometry.Ny-1 ; j>=0; j--)
                {
                    if(i==geometry.Nx-1)
                    {
                        psi_right_j[j] = Storage_R_boun_psi[g][j][n][0];
                    }
                    else psi_right_j[j] = psi_left_j[j];

                    if(j==geometry.Ny-1)
                    {
                        psi_up_i[i] = Storage_T_boun_psi[g][i][n][0];
                    }
                    else psi_up_i[i] = psi_down_i[i];

                    psi_mid = 1.0/(geometry.sigma_t(i,j, g)+ 2.0 * mu /geometry.delx_i[i] + 2.0 * eta/geometry.dely_j[j])
                              *(+2.0*mu/geometry.delx_i[i] * psi_right_j[j] + 2.0*eta/geometry.dely_j[j]*psi_up_i[i]+ q[i][j]) ;

                    phi[i][j] = phi[i][j]+0.25 * w*psi_mid;

                    psi_left_j[j] =  (2.0 * psi_mid - psi_right_j[j]);
                    psi_down_i[i] =  (2.0 * psi_mid - psi_up_i[i]);

                    if(i==0)
                    {
                        Storage_L_boun_psi[g][j][n][1] = geometry.alpha_L*psi_left_j[j];
                    }
                    if(j==0)
                    {
                        Storage_B_boun_psi[g][i][n][1]= geometry.alpha_B*psi_down_i[i];
                    }
                }
            }

        }

        for(int n=0; n<angle.total_num; n++)//mu>0 and eta<0 + -
        {
            mu = angle.mu[angle.mu_sequence[n]];
            eta = angle.mu[angle.eta_sequence[n]];
            w = angle.w[angle.sequence[n]];

            for(int i=0; i<geometry.Nx; i++)
            {
                for(int j=geometry.Ny-1; j>=0; j--)
                {
                    if(i==0)
                    {
                        psi_left_j[j] = Storage_L_boun_psi[g][j][n][1];
                    }
                    else psi_left_j[j]=psi_right_j[j];

                    if(j==(geometry.Ny-1))
                    {
                        psi_up_i[i] = Storage_T_boun_psi[g][i][n][1];
                    }
                    else psi_up_i[i] =psi_down_i[i];

                    psi_mid = 1.0/(geometry.sigma_t(i,j, g)+ 2.0 * mu /geometry.delx_i[i] + 2.0 * eta/geometry.dely_j[j])
                              *(2.0*mu/geometry.delx_i[i] * psi_left_j[j] + 2.0*eta/geometry.dely_j[j]*psi_up_i[i]+ q[i][j]) ;


                    psi_right_j[j] =  (2.0 * psi_mid - psi_left_j[j]);
                    psi_down_i[i] =(2.0 * psi_mid - psi_up_i[i]);

                    phi[i][j] = phi[i][j]+0.25 * angle.w[angle.sequence[n]]*psi_mid;

                    if(i==geometry.Nx-1)
                    {
                        Storage_R_boun_psi[g][j][n][0] = geometry.alpha_R*psi_right_j[j];
                    }
                    if(j==0)
                    {
                        Storage_B_boun_psi[g][i][n][0]= geometry.alpha_B*psi_down_i[i];
                    }
                }
            }
        }



        return 1;
    };

};





double sum_fission_rate(vector<vector<vector<double>>> const &flux_g_ij, geometry_class const & geometry)
{
    double fission_rate=0.0;
    for (int i=0; i<geometry.Nx; i++)
    {
        for(int j=0; j< geometry.Ny; j++)
        {
            for (int g=0; g<geometry.groups; g++)
                fission_rate = fission_rate + flux_g_ij[g][i][j] * geometry.nu_sigmaf(i,j,g) * geometry.delx_i[i] * geometry.dely_j[j];

        }
    }

    return fission_rate;

}



int main()
{

    
    int mesh_num;
    cout << "give number of mesh in one axis"<<endl;
    cin>>mesh_num;
    //angular discretization
    angular angle;
    //form geometry
    geometry_class geometry(mesh_num, r_c);

    double keff = 1.0, keff_old=1.0;
    Transport_sweep_class sweep_object(angle, geometry);

    vector<vector<vector<double>>> flux_g_ij(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));
    vector<vector<vector<double>>> flux_g_in(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));
    vector<vector<vector<double>>> flux_g_old(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));


    vector<vector<double>> q (geometry.Nx, vector<double> (geometry.Ny,1.0));
    vector<vector<vector<double>>> fission_density_g_ij(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));


    cout<< "initialization done"<<endl;

    ofstream outputfile;
    string str = name;
    str+= ".out";
    outputfile.open(str, ios::app);

    int out_it = 0;
    int src_it = 0;
    int total_src_it = 0;
    double scatter_density = 0;
    outputfile<<"iteration \t keff \t \t error\t"<<endl;

    do //outer iteration
    {
        geometry.calculate_fission_density_g(flux_g_ij, fission_density_g_ij, keff);
        keff_old = keff;
        cout<<"here at outer iteration"<<endl;
    
        src_it = 0;
        do//inner
        {
            vector3Dcopy(flux_g_in, flux_g_ij);
            for(int g=0; g<geometry.groups; g++)
            {
                for(int i=0; i<geometry.Nx; i++)
                {
                    for (int j=0; j<geometry.Ny; j++)
                    {
                        scatter_density = 0;
                        for(int from_g=0; from_g<geometry.groups; from_g++)
                        {
                            scatter_density = scatter_density + flux_g_ij[from_g][i][j] * geometry.sigma_s(i,j,g,from_g);       
                        }
                        q[i][j] = scatter_density + fission_density_g_ij[g][i][j];
                    }

                }
                //cout<<abs_max(fraction_error(flux_g_ij[g], flux_g_old[g]))<<endl;
                
                sweep_object.do_sweep(geometry, angle, q, g, flux_g_ij[g]);
                
            } 
            src_it+=1;
        } while (!is_flux_converged(flux_g_in, flux_g_ij) && src_it<= max_it);
    
        total_src_it = total_src_it + src_it;

        keff = keff_old * sum_fission_rate(flux_g_ij, geometry)/sum_fission_rate(flux_g_old, geometry);
        vector3Dcopy(flux_g_old, flux_g_ij);
        out_it +=1;
        outputfile<<out_it<<"\t \t"<<keff<<"\t \t"<<keff-keff_old<<endl;
        cout<<out_it<<"\t"<<keff<<"\t \t"<<keff-keff_old<<endl;

    }while(!is_keff_converged(keff_old, keff));

    cout<<"outer iteration done: "<<keff << " outer iteration: "<< out_it <<" source iteration: " << total_src_it << endl;
    cout<<endl<<endl;
    


    // if (outputfile.is_open()) {
    //     outputfile<<name<<"  delx = "<< geometry.delx_i[0] << " keff= " << keff << " flux_ratio: "<< flux_g_ij[0][0][0]/flux_g_ij[1][0][0];
    //     outputfile<<endl<<endl;
    //     for(int g= 0 ;  g<geometry.groups; g++)
    //     {
    //         outputfile << "group "<< g<<endl<<endl;

    //         for (int j = geometry.Ny-1; j>=0; j--)
    //         {

    //             for(int i = 0; i<geometry.Nx; i++)
    //             {
    //                 outputfile<<flux_g_ij[g][i][j]<< " ";
    //                 if(i==geometry.Nx)
    //                     outputfile<< ";"<<endl;
    //             }              
    //         }
    //         outputfile << endl;
    //     }

    //     outputfile << endl <<endl;

        

    // }

    if (outputfile.is_open()) 
    {
        outputfile<<name<<"  delx = "<< geometry.delx_i[0] << " keff= " << keff << " with error: "<< keff-keff_old;
        outputfile<<endl<<endl;
        outputfile<<"pos"<<"\t flux1 \t \t flux2"<<endl;
        for(int i=0;i<geometry.Nx; i++)
        {
            outputfile<<i<<"\t";
            for(int g=0; g<geometry.groups; g++)
            {
                outputfile<< flux_g_ij[g][i][0]<< "\t \t";
            }
            outputfile<<endl;
        }

        outputfile << endl <<endl;

    }

    else {
            std::cerr << "Unable to open file!" << std::endl;
            return 1;
        }

    return 0;

}