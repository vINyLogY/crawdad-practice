#include <iostream>
#include <string>
#include <cstdio>

#include "atomdata.h"
#include "SCF.h"
//#include "Molecule.h"

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

    Basis basis("enuc.dat", "s.dat", "t.dat", "v.dat", "eri.dat"); 

    basis.iteration(5, 1e-5);

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
        ij = INDEX((i-1), (j-1));
        kl = INDEX((k-1), (l-1));
        ijkl = INDEX(ij, kl);
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

void Basis::set_fx()
{
    if(iter == 0) F = H;
    else
    {
        int ij, kl, ijkl, ik, jl, ikjl;
        F = H;
        for(int i=0; i < max; i++)
            for(int j=0; j < max; j++)
            {
                for(int k=0; k < max; k++)
                    for(int l=0; l < max; l++)
                    {
                        ij = INDEX(i, j);
                        kl = INDEX(k, l);
                        ijkl = INDEX(ij, kl);
                        ik = INDEX(i, k);
                        jl = INDEX(j, l);
                        ikjl = INDEX(ik, jl);
                        F(i, j) += D(k, l) * (2.0 * tei[ijkl] - tei[ikjl]);
                     }
            }
    }
    Fx = Sx.adjoint() * F * Sx;

    return;
}

void Basis::set_eps()
{
    Eigen::SelfAdjointEigenSolver<Matrix> solver(Fx);
    Cx = solver.eigenvectors();
    Matrix e = solver.eigenvalues();
    Eps.resize(max, max);
    for(int j = 0 ; j < max; j ++)
    {
        Eps(j, j) = e(j);
    }

    return;
}

void Basis::set_d()
{
    if(iter != 0)
    {
        _D = D;
        _Eele = Eele;
    }

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
    Eele = 0.0;
    for(int i=0; i < max; i++)
        for(int j=0; j < max; j++)
        {
            Eele += D(i, j) * (H(i, j) + F(i, j));
        }
    Etot = Eele + Etot;
    return;
}

bool Basis::e_conv(double delta)
{
    if(Eele - _Eele < delta && Eele - _Eele > 0.0 - delta) 
        return true;
    else return false;

}

bool Basis::rms_conv(double delta)
{
    double ans = 0;
    for(int i = 0; i < max; i++)
        for(int j = 0; j < max; j++)
            ans += (D(i, j) - _D(i, j)) * (D(i, j) - _D(i, j));
    ans = sqrt(ans);

    if(ans < delta) 
        return true;
    else return false;        
}

void Basis::iteration(int oao, double err)
{
    occ = oao;
    iter = 0;

    while(true)
    {
        set_fx();
        set_eps();
        set_d();
        printf("Iter %2d E(elec) = %12.8lf\n", iter, Eele);
        iter += 1;
        if(e_conv(err) && rms_conv(err)) break;
    }

    return;
}


