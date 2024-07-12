#include <iostream>
#include<vector>

using namespace std;
int foo(vector<vector<double>> &psi_ij)
{
    psi_ij[1][2] = 68;
    cout<<psi_ij[1][2]<<endl;
    return 0;
}

int main()
{
    int Nx = 10, Ny=5;
    vector<vector<double>> psi_ij(Nx,vector<double> (Ny,0));
    vector<vector<double>> psi_ij_half(Nx+1, vector<double> (Ny+1,0) );
    for (int i=0; i<size(psi_ij);i++)
    {
        for (int j=0; j<size(psi_ij[0]);j++)
        {

            cout<<psi_ij[i][j]<<"\t";
        }
        cout<<endl;
    }
    foo(psi_ij);

    for (int i=0; i<size(psi_ij);i++)
    {
        for (int j=0; j<size(psi_ij[0]);j++)
        {
            
            cout<<psi_ij[i][j]<<"\t";
        }
        cout<<endl;
    }
    
}