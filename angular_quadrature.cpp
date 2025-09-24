#include "angular_quadrature.hpp"

angular::angular(int N)
{
    if(N==2) //S2 quadrature
    {
        //mapping weights for each angluar direction
        sequence = std::vector<int> { 
            1
        };

        mu = std::vector<double> {
            1.0/sqrt(3.0)
        };

        w = std::vector<double> {1.0};
    }

    if(N==4) //S4
    {
        //mapping weights for each angluar direction
        sequence = std::vector<int> {
            1,
            1, 1
        };

        mu = std::vector<double> {
            0.30163878,
            0.90444905
        };

        w = std::vector<double> {0.33333333};
    }

    if(N==6) // S6
    {
        //mapping weights for each angluar direction
        sequence = std::vector<int>{
            1,
            2, 2,
            1, 2, 1
        };
        mu = std::vector<double>{
            0.23009194,
            0.68813432,
            0.94557676
        };
        w = std::vector<double> {
            0.16944656,
            0.16388677
        };
    }

    if(N==8) //S8
    {
        //mapping weights for each angluar direction
        sequence = std::vector<int> {
            1,
            2, 2,
            2, 3, 2,
            1, 2, 2, 1
        };
        mu = std::vector<double> { 
            0.19232747,
            0.57735027,
            0.79352178,
            0.96229948
        };
        w = std::vector<double> 
        {
            0.11678847,
            0.09325523,
            0.09010320
        };
    }
    
    if(N==12) //S12
    {
        //mapping weights for each angluar direction
        sequence = std::vector<int> {
            1,
            2, 2,
            3, 4, 3,
            3, 5, 5, 3,
            2, 4, 5, 4, 2,
            1, 2, 3, 3, 2, 1
        };
        mu=std::vector<double>{
            0.15395746,
            0.45769112,
            0.62869660,
            0.76225828,
            0.87568027,
            0.97600932
        };
        w = std::vector<double>{
            0.07332178,
            0.05266740,
            0.04161495,
            0.03895667,
            0.03249018
        };
    }

    if(N==16) //S16
    {
        std::vector<int> sequence1{ 
            1,
            2, 2,
            3, 5, 3,
            4, 6, 6, 4,
            4, 7, 8, 7, 4,
            3, 6, 8, 8, 6, 3,
            2, 5, 6, 7, 6, 5, 2,
            1, 2, 3, 4, 4, 3, 2, 1 
        };

        sequence = sequence1;
        std::vector<double>mu1{  0.13344572,
                                0.39119433,
                                0.53689687,
                                0.65075610,
                                0.74746822,
                                0.83302700,
                                0.91058181,
                                0.98203079};
        mu = mu1;
        std::vector<double>w1{  0.05415425,
                                0.03679653,
                                0.02777273,
                                0.02580284,
                                0.02494275,
                                0.01962325,
                                0.01879762,
                                0.01544801 };
        w=w1;
    }
    
    //mapping mu and eta for each angular direction
    for(int i=0; i<std::size(mu) ; i++)
    {
        for(int j=0;j<=i;j++)
        {
            mu_sequence.push_back(i-j);
            eta_sequence.push_back(j);
        }
    }

    // This loop subtracts 1 from all sequence values to make them 0-indexed
    for (int i=0; i<std::size(sequence);i++)
        sequence[i]=sequence[i]-1;

    total_num = std::size(sequence);
}