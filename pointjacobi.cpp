#include <iostream>
#include <cmath>
#include <fstream>

float u[41][31], u_[41][31] ; // storing value of grid where, u= psi
float plot[41][31];           // final grid value stored used for plotting 
const int value = 100; // 100, 150, 200 are the given initial streamlines values
float e = 0;
float k = 0;
float error = 0.0;
int i, j;
float norm, norm_difference; // used to calculate error 
using namespace std;

class pointjacobi
{
public:
    void boundarycondn();
    void initialcondn();
    void function();
};

void pointjacobi::initialcondn()  // storing initial values in the matrix initially
{
    for (int t = 1; t < 40; t++)   // index of rows of the grid 
    {
        for (int n = 1; n < 30; n++)  // index of column of the grid
        {

            u[t][n] = value;

            // cout << u[t][n] << " ";
        }
        // cout << endl;
    }
}

void pointjacobi::boundarycondn()   // given boundary condition 
{
    for (int t = 0; t < 41; t++) // index of rows of the grid 
    {
        for (int n = 0; n < 31; n++)  /// index of column of the grid
        {
            if (t == 40 && n < 12)
            {
                u[t][n] = 100;
                u_[t][n] = 100;
            }
            else if (t == 40 && n > 11)
            {
                u[t][n] = 200;
                u_[t][n] = 200;
            }
            else if (t == 0)
            {
                u[t][n] = 100;
                u_[t][n] = 100;
            }
            u[t][0] = 100;
            u_[t][0] = 100;
            if (t < 22)
            {
                u[t][30] = 100;
                u_[t][30] = 100;
            }
            else
            {
                u[t][30] = 200;
                u_[t][30] = 200;
            }
            // cout << u[t][n] << " ";
        }
        cout << endl;
    }
}

void pointjacobi::function()  // using point jacobi we calculate the values at each grid 3*4
{
    for (int i = 0; i < 1000; i++)    // i is the number of the iterations carried out
    {
        for (int t = 1; t < 40; t++)
        {
            for (int n = 1; n < 30; n++)
            {
                u_[t][n] = (u[t - 1][n] + u[t][n + 1] + u[t][n - 1] + u[t + 1][n]) * 0.25;
            }
        }

        e = 0.0;
        k = 0.0;

        for (int t = 1; t < 40; t++)
        {
            for (int n = 1; n < 30; n++)
            {
                e = pow(u[t][n], 2) + e;                      // norm
                k = pow(abs(u_[t][n]) - abs(u[t][n]), 2) + k; // norm difference
                u[t][n] = u_[t][n];
            }
        }

        error = pow((k / e), 0.5);    // error calculation

          // condition as given in question
        if (error < 0.0001)
        {
            cout << "error: " << error << " " << endl;
            cout << "Iterations: " << i << endl;

            // printing the final matrix
            for (int t = 0; t < 41; t++)
            {
                for (int n = 0; n < 31; n++)
                {
                    cout << u[t][n] << " ";
                }
                cout << endl << endl;
            }
              // saving the converged matrix 
            ofstream outfile("plot100.txt");
            for (int t = 0; t < 41; t++)
            {
                for (int n= 0; n< 31; n++)
                {
                    outfile << u_[t][n] << " ";
                }
                outfile << "\n";
            }
            outfile.close();

            break;
        }
    }
}
int main()
{
    pointjacobi pj;
    pj.initialcondn();
    pj.boundarycondn();
    cout << "Printing Final Matrix:" << endl;
    pj.function();
    return 0;
}