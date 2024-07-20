#include <iostream>
#include<vector>
#include<cmath>

using namespace std;

int main()
{
    vector<vector<vector<vector<double>>>> bottom_top(2, vector<vector<vector<double>>> (4, vector<vector<double>> (10, vector<double> (3,13))));
    cout <<(bottom_top[1][2][2][2]);
    string str="bello";
    cout<<str;

    vector<vector<double>> flux_old(10, vector<double> (20,0.0));
    vector<vector<double>> q(10, vector<double> (20,0.0));
    vector<vector<double>> fission_source = q, fission_source_new;
    for(int i=0; i<10; i++)
    {
        for (int j=0; j<20; j++)
        {
            flux_old[i][j] = 1.0;
            fission_source[i][j] = flux_old[i][j] * 0.86;
        }
        
    }
    
    if( 0.1 < 10^(-5) )
    {
        double a = pow(10,-5);
        cout<<a;
    }
    
}