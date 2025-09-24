#ifndef ANGUALR_QUADRATURE
#define ANGULAR_QUADRATURE
#include<vector>
#include<cmath>

class angular
{
    public:
    std::vector<int> sequence;
    std::vector<double> mu, mu_sequence, eta_sequence;
    std::vector<double> w;
    int total_num;
    angular(int N);
    
};

#endif