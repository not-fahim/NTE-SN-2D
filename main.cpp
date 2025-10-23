#include <iostream>
#include<vector>
#include <cctype>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<sstream>
#include "geometry_class.hpp"
#include "angular_quadrature.hpp"
#include "geo_mat_parser.hpp"

using namespace std;

double intol, outtol;
int max_it;

void print_banner(ofstream &outputfile) {
    // Copy and paste the text from one of the banners below into this string
    std::string banner_text = R"(
**************************************************
*         .-~*~~~*~-.                           *
*        .-~~~~~~~~~-.                           *
*       /  X       X  \                          *
*      |    .-----.    |                         *
*       \  '_______'  /                          *
*        `-.........-'                           *
* *  ____   ___  ____ ____                       *
*   |  _ \ / _ \_   _| __ )                      *
*   | | | | | | || | |  _ \                      *
*   | |_| | |_| || | | |_) |                     *
*   |____/ \___/ |_| |____/                      *
* *
*   Discrete Ordinate with Two Braincells        *
**************************************************
    )";
    outputfile << banner_text << std::endl;
}

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


double Linferror(vector<vector<vector<double>>> const &vector1, vector<vector<vector<double>>> const &vector2)
{

        double abs_max=0;
        double a;

        for (int i=0; i< size(vector1); i++)
        {
            for (int j=0; j< size(vector1[i]); j++)
            {
                for(int k=0; k<size(vector1[i][j]); k++)
                {
                    if(i==0 && j==0 && k==0) 
                        abs_max = abs(vector1[i][j][k] - vector2[i][j][k]);
                
                    a=abs(vector1[i][j][k] - vector2[i][j][k]);

                    if(a>abs_max) 
                        abs_max = a; //inf norm
                } 
            }
        }

        return abs_max;

}

template <typename t>
bool is_flux_converged(vector<vector<vector<t>>> const &vector1, vector<vector<vector<t>>> const &vector2 )
{
    if(Linferror(vector1, vector2) < intol )
    {
        return true;
    }
    else
    {
        return false;
    }
}


double Linferror(vector<vector<double>> const &vector1, vector<vector<double>> const &vector2)
{
        double abs_max=0;
        double a;


        for (int i=0; i< size(vector1); i++)
        {
            for (int j=0; j< size(vector1[0]); j++)
            {
                if(i==0 && j==0)
                    abs_max = (vector1[i][j] - vector2[i][j]);
                
                a=abs(vector1[i][j] - vector2[i][j]);

                if(a>abs_max) 
                    abs_max = a; 
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
bool is_flux_converged(vector<vector<t>> const &vector1, vector<vector<t>> const &vector2 )
{
    if(Linferror(vector1, vector2) < intol )
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

    input_class input_object = read_input_file();
    print_input_data(input_object);
    intol = input_object.tol_in;
    outtol = input_object.tol_out;
    max_it = input_object.max_it;

    //form angle
    angular angle(input_object.S_n);
    angle.print_quadrature_card();
    //form geometr_n
    geometry_class geometry(input_object);

    double keff = 1.0, keff_old=1.0;
    Transport_sweep_class sweep_object(angle, geometry);

    vector<vector<vector<double>>> flux_g_ij(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));
    vector<vector<vector<double>>> flux_g_in(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));
    vector<vector<vector<double>>> flux_g_old(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));


    vector<vector<double>> q (geometry.Nx, vector<double> (geometry.Ny,1.0));
    vector<vector<vector<double>>> fission_density_g_ij(geometry.groups, vector<vector<double>> (geometry.Nx, vector<double> (geometry.Ny,1.0)));


    cout<< "initialization done"<<endl;

    ofstream outputfile;
    string name = input_object.name;
    string str = input_object.name;
    str+= ".out";
    outputfile.open(str, ios::app);

    int out_it = 0;
    int src_it = 0;
    int total_src_it = 0;
    double scatter_density = 0;
    print_banner(outputfile);
    outputfile << std::fixed << std::setprecision(8);
    outputfile<<name<<"\t"<<geometry.groups<<" group"<< endl;
    outputfile<<"XS file: "<<input_object.xsfile<<endl;
    outputfile<<"S-"<<input_object.S_n<< "\t quad set, \t mesh refinement: "<<input_object.refinement<<endl <<endl;
    outputfile<<"iteration \t keff \t \t error \t \t iteration "<<endl;
    cout<<"iteration \t keff \t \t error \t iteration"<<endl;

    do //outer iteration
    {
        geometry.calculate_fission_density_g(flux_g_ij, fission_density_g_ij, keff);
        keff_old = keff;

    
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
        outputfile<<out_it<<"\t \t"<<keff<<"\t \t"<<(keff-keff_old)/keff<< "\t \t"<< src_it<<endl;
        cout<<out_it<<"\t"<<keff<<"\t \t"<<(keff-keff_old)/keff<<"\t \t"<< src_it <<endl;

    }while(!is_keff_converged(keff_old, keff) );

    cout<<"\n ---outer iteration done---\n keff: "<<keff << " outer iteration: "<< out_it <<" source iteration: " << total_src_it << endl;
    cout<<endl<<endl;
    

    if (outputfile.is_open()) 
    {
        outputfile<<name<<"  refinement = "<< input_object.refinement<<" S-"<<input_object.S_n << " keff= " << keff ;
        outputfile<<endl;
        outputfile<< "flux distribution is saved in the flux_"<<name<<".m and .py file script"<< endl;
        outputfile<<"--- End of calculation ---"<<endl<<endl;
    }

    else 
    {
        std::cerr << "Unable to open file!" << std::endl;
        return 1;
    }
    outputfile.close();
    return 0;

}