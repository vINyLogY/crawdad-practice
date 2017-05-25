#include "Molecule.h"
#include "atomdata.h"

#include <iostream>
#include <cstdio>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;
 
//截断长度
#define LEN 4.0

#define DEBUG

using namespace std;

int main(void)
{
    #ifdef DEBUG
        freopen("input.txt", "r", stdin);
        freopen("output.txt", "w", stdout);
    #endif

    Molecule mol("project1.dat");

    printf("Number of atoms: %d\n", mol.natom);

    printf("\nInput Cartesian coordinates:\n");
    mol.print_geom();

    printf("\nInteratomic distances (bohr):\n");
    for(int i = 1; i < mol.natom; i++)
        for(int j = 0; j < i; j++)
            printf("%d %d %8.5f\n", i, j, mol.bond(i, j));

    printf("\nBond angles:\n");
    for(int j = 0; j < mol.natom; j++)
        for(int i = 0; i < mol.natom; i++)
            for(int k = 0; k < i; k++)
                if(i != j && k != j
                    &&mol.bond(i, j) < LEN && mol.bond(j, k) < LEN)
                    printf("%d-%d-%d %12.6f\n", i, j, k, mol.angle(i, j, k) * 180.0 / PI);
                        
    printf("\nOut-of-plane angles:\n");
    for(int l = 0; l < mol.natom; l++)
        for(int j = 0; j < mol.natom; j++)
            for(int i = 0; i < mol.natom; i++)
                for(int k = 0; k < i; k++)
                    if(l != i && l != j && l != k && i != j && j != k
                        && mol.bond(l, j) < LEN && mol.bond(i, j) < LEN && mol.bond(k, j) < LEN)
                        printf("%d-%d-%d-%d %12.6f\n", l, i, j, k, mol.angle_op(l, i, j, k) * 180.0 / PI);

    printf("\nTorsional angles:\n");
    for(int j=0; j < mol.natom; j++)
        for(int k=0; k < j; k++)
            for(int i=0; i < mol.natom; i++)
                for(int l=0; l < mol.natom; l++)
                    if(i != j && i != k && l != j && l != k && i != l
                        && mol.bond(i,j) < LEN && mol.bond(j,k) < LEN && mol.bond(k,l) < LEN)
                        printf("%d-%d-%d-%d %12.6f\n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/PI));

    printf("\nMolecular center of mass:");
    for(int i = 0; i < DIM; i++)
        printf(" %12.6f", mol.com(i));
    printf("\n");
    mol.translate(0.0 - mol.com(0), 0.0 -mol.com(1), 0.0 -mol.com(2));

    Matrix I(3,3);
    double mi;

    for(int i=0; i < mol.natom; i++)
    {
        mi = masses[mol.zvals[i]];
        I(0,0) += mi * (mol.geom[i][1]*mol.geom[i][1] + mol.geom[i][2]*mol.geom[i][2]);
        I(1,1) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][2]*mol.geom[i][2]);
        I(2,2) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][1]*mol.geom[i][1]);
        I(0,1) += mi * mol.geom[i][0]*mol.geom[i][1];
        I(0,2) += mi * mol.geom[i][0]*mol.geom[i][2];
        I(1,2) += mi * mol.geom[i][1]*mol.geom[i][2];
    }
    
    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);
    
    cout << "\nMoment of inertia tensor (amu bohr^2):\n";
    cout << I << endl;
    
    // find the principal moments
    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();
    
    cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
    cout << evals << endl;
    
    double conv = 0.529177249 * 0.529177249;
    cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
    cout << evals * conv << endl;
    
    conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
    cout << "\nPrincipal moments of inertia (g * cm^2):\n";
    cout << evals * conv << endl;
    
    // classify the rotor 
    if(mol.natom == 2) cout << "\nMolecule is diatomic.\n";
    else if(evals(0) < 1e-4) cout << "\nMolecule is linear.\n";
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
        cout << "\nMolecule is a spherical top.\n";
    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) > 1e-4))
        cout << "\nMolecule is an oblate symmetric top.\n";
    else if((fabs(evals(0) - evals(1)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
        cout << "\nMolecule is a prolate symmetric top.\n";
    else cout << "\nMolecule is an asymmetric top.\n";

    // compute the rotational constants 
    const double _pi = acos(-1.0);
    conv = 6.6260755E-34/(8.0 * _pi * _pi);
    conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
    conv *= 1e-6;
    cout << "\nRotational constants (MHz):\n";
    cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;
    
    conv = 6.6260755E-34/(8.0 * _pi * _pi);
    conv /= 1.6605402E-27 * 0.529177249E-10 * 0.529177249E-10;
    conv /= 2.99792458E10;
    cout << "\nRotational constants (cm-1):\n";
    cout << "\tA = " << conv/evals(0) << "\t B = " << conv/evals(1) << "\t C = " << conv/evals(2) << endl;
    
    return 0;
}
