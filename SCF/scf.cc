#include <iostream>
#include <string>
#include <cstdio>

#include "atomdata.h"
#include "SCF.h"
#include "Molecule.h"

#include "Eigen/Dense"

//#define DEBUG

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

//STEP 1: Nuclear Repulsion Energy
//STEP 2: One-Electron Integrals
//STEP 3:Two-Electron Integrals

    Basis basis("enuc.dat", "s.dat", "t.dat", "v.dat", "eri.dat");

//STEP 4: Build the Orthogonalization Matrix

//    cout << basis.Sx << endl;

//STEP 5: Build the Initial (Guess) Density

    basis.set_fx(0);
    basis.set_eps(0);
    basis.occ = 5;
    basis.set_d(0);
    cout << basis.D << endl;
    basis.Eele = 0.0;
    for(int i = 0; i < basis.occ; i++)
    {
        basis.Eele += 2 * basis.Eps(i, i);
    }
    basis.Etot = basis.Eele + basis.Enuc;
    cout << basis.Etot <<endl;




//STEP 6: Compute the Initial SCF Energy

//STEP 7: Compute the New Fock Matrix

//STEP 8: Build the New Density Matrix

//STEP 9: Compute the New SCF Energy

//STEP 10: Test for Convergence

    return 0;
}


Basis::Basis(const char *fenuc, const char *fs, const char *ft, const char *fv, const char *feri)
{
    //read the nuclear repulsion energy.
    FILE *f0 = fopen(fenuc, "r");
    fscanf(f0, "%lf", &Enuc);
    fclose(f0);

    //read the one-electron integrals.
    S = fr_basis(fs);
    T = fr_basis(ft);
    V = fr_basis(fv);
    H = T + V;
    max = S.cols();

    //read the two-electron integrals.
    tei = new double[max * max * max * max];
    FILE *f2 = fopen(feri, "r");
    int i, j, k, l, ij, kl, ijkl;
    double temp;
    while(fscanf(f2, "%d %d %d %d %lf", &i, &j, &k, &l, &temp) != EOF)
    {
        ij = (i > j) ? (i*(i+1)/2)+j : (j*(j+1)/2+i);
        kl = (k > l) ? (k*(k+1)/2)+l : (l*(l+1)/2+k);
        ijkl = (ij > kl) ? (ij*(ij+1)/2+kl) : (kl*(kl+1)/2+ij);
        tei[ijkl] = temp; 
    }
    fclose(f2);

    //Orthogonalization
    Eigen::SelfAdjointEigenSolver<Matrix> solver(S);
    Matrix Ls = solver.eigenvectors();
    Matrix lambda = solver.eigenvalues();
    Matrix La;
    La.resize(max, max);
    for (int i = 0; i < max; i++)
        La(i, i) = 1 / sqrt(lambda(i));
    Sx = Ls * La * Ls.adjoint();
}

Basis::~Basis()
{
    delete[] tei;
}

Matrix fr_basis(const char *filename)
{
    Matrix ans;
    double temp;
    int i, j, max;
    FILE *f = fopen(filename, "r");
    while(fscanf(f, "%d %d %lf", &i, &j, &temp) != EOF);
    max = i;
    rewind(f);
    ans.resize(max, max);
    while(fscanf(f, "%d %d %lf", &i, &j, &temp) != EOF)
    {
        ans(i - 1, j - 1) = temp;
    }
    for(int _i = 0; _i < max - 1; _i++)
        for(int _j = _i + 1; _j < max; _j++)
            ans(_i, _j) = ans(_j, _i);
    return ans;
}

void Basis::set_fx(int i)
{
    if(i == 0)
        Fx = Sx.adjoint() * H * Sx;
}

void Basis::set_eps(int i)
{
    if(i == 0)
    {
        Eigen::SelfAdjointEigenSolver<Matrix> solver(Fx);
        Cx = solver.eigenvectors();
        Matrix e = solver.eigenvalues();
        Eps.resize(max, max);
        for(int j = 0 ; j < max; j ++)
        {
            Eps(j, j) = e(j);
        }
        cout << "e = " << e <<endl;
    }
}

void Basis::set_d(int i)
{
    if(i == 0)
    {
        C = Sx * Cx;
        D.resize(max, max);
        double temp;
        for(int i = 0; i < max; i++)
            for(int j = 0; j < max; j++)
            {
                temp = 0.0;
                for(int m = 0; m < occ; m++)
                    temp += C(i,m) * C(j,m); 
                D(i,j) = temp;
            }
    }
}
