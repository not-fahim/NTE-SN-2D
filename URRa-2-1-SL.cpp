char name[]="URRa-2-1-SL";

double sigma_s1_g[2][2] = {{0.27459, 0.0},
                           {0.83318, 0.0075737}};
double sigma_t_fuel[2] = {0.65696-sigma_s1_g[0][0]/3, 2.52025-1/3*(sigma_s1_g[1][0]+sigma_s1_g[1][1])};

double sigma_s_fuel[2][2] = {{0.62568-1/3*sigma_s1_g[0][0], 0.0},
                                {0.029227, 2.44383-1/3*(sigma_s1_g[1][0]+sigma_s1_g[1][1])}
                                };

double nu_sigmaf_fuel[2] = {2.50*0.0010484, 2.50*0.050632 };
double chi_fuel[2] = {1, 0};

double BC[]={1,1,1,0}; //B T L R
double r_c = 9.4959;