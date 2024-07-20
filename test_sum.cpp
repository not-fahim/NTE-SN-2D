#include<iostream>
#include<vector>

using namespace std;
template <typename t>
t sum(vector<vector<t>> vector1)
{
    t sum;
    for (int i=0; i<size(vector1); i++)
    {
        for(int j=0; j<size(vector1[0]); j++)
        {
            sum = sum+vector1[i][j];
        }
    }

    return sum;
}

int main()
{
    vector<vector<double>> a ={{3,4,5}, {4,5,6}, {2,3,4}};
    cout<<sum(a);
}