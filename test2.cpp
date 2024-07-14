#include <iostream>
#include<vector>

using namespace std;

int main()
{
    vector<vector<vector<vector<double>>>> bottom_top(2, vector<vector<vector<double>>> (4, vector<vector<double>> (10, vector<double> (3,13))));
    cout <<(bottom_top[1][2][2][2]);
    string str="bello";
    cout<<str;
    
}