#include "Molecule.h"
#include "atomdata.h"

#include <iostream>
#include <cstdio>
#include <cmath>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
 
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
 
//截断长度
#define LEN 4.0

#define DEBUG

using namespace std;

int  main(void)
{
    #ifdef DEBUG
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
    #endif

    Molecule mol("project2.dat");
    mol.read_mw_hessian("hessian.dat");
    Eigen::SelfAdjointEigenSolver<Matrix> solver(mol.H);
    Matrix lambda = solver.eigenvalues();
    const double con = 5140.485;
    for(int i = mol.natom * 3 - 1; i >= 0; i--)
    {
        if(lambda(i) < 0)
            lambda(i) = 0;
        printf("%d %12.6f\n", i, con * sqrt(lambda(i)));
    }

    return 0;
}