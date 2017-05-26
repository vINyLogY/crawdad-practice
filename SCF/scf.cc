#include <iostream>
#include <string>
#include <cstdio>

#include "atomdata.h"
#include "SCF.h"
//#include "Molecule.h"

#include "Eigen/Dense"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

int main(void)
{
    #ifdef DEBUG
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
    #endif

//    Molecule mol("geom.dat");
    double enuc = fr_Enuc("enuc.dat");

/*
    Matrix S = fr_basis("s.dat"),
           T = fr_basis("t.dat"),
           V = fr_basis("v.dat"),
           H;
    H = T + V;
    cout << T << endl;
    cout << V << endl;
    cout << H << endl;
*/
    Matrix Ans;
    int i, j, max;
    double temp;
    FILE *f = fopen("s.dat", "r");
    while(fscanf(f, "%d %d %lf", &i, &j, &temp) != EOF);
    max = i;
    rewind(f);
    Ans.resize(max, max);
    while(fscanf(f, "%d %d %lf", &i, &j, &temp) != EOF)
    {
        Ans(i - 1, j - 1) = temp;
    }
    for(int _i = 0; _i < 6; _i++)
        for(int _j = _i + 1; _j < 7; _j++)
        {
            Ans(_i, _j) = Ans(_j, _i);
        }
    
    cout << Ans << endl;

    fclose(f);


    


    return 0;
}

double fr_Enuc(const char *filename)
{
    double enuc;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%lf", &enuc);
    fclose(f);
    return enuc;
}
