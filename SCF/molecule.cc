#include "Molecule.h"
#include "atomdata.h"

#include <iostream>
#include <cstdio>
#include <cmath>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

//Project 1
void Molecule::print_geom()
{
    for(int i = 0; i < natom; i++)
        printf("%d %10.6f %10.6f %10.6f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::rotate_z(double phi)
{
    double x, y;
    for(int i = 0; i < natom; i++)
    {
        x = geom[i][0];
        y = geom[i][1];
        geom[i][0] = cos(phi * PI / 180.0) * x - sin(phi * PI / 180.0) * y;
        geom[i][1] = sin(phi * PI / 180.0) * x + cos(phi * PI / 180.0) * y;

    }
}
 
void Molecule::translate(double x, double y, double z)
{
    for(int i = 0; i < natom; i++)
    {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

double Molecule::bond(int atom1, int atom2)
{
    double ans, v1[DIM];
    for(int i = 0; i < DIM; i++)
        v1[i] = geom[atom1][i] - geom[atom2][i];
    ans = sqrt(dot(v1, v1));
    return ans;

}

double Molecule::angle(int a, int b, int c)
{
    double ans, v1[DIM], v2[DIM];
    for(int i = 0; i < DIM; i++)
    {
        v1[i] = geom[a][i] - geom[b][i];
        v2[i] = geom[c][i] - geom[b][i]; 
    }
    ans = dot(v1, v2)/bond(a, b)/bond(c, b);
    if(ans < -1.0) ans = acos(-1.0);
    else if(ans > 1.0) ans = acos(1.0);
    else ans = acos(ans);
    return ans;
    
}

double Molecule::angle_op(int d, int a, int b, int c)
{
    double v0[DIM], v1[DIM], v2[DIM], v3[DIM], ans;
    for(int i = 0; i < DIM; i++)
    {
        v0[i] = geom[d][i] - geom[b][i];
        v1[i] = geom[a][i] - geom[b][i];
        v2[i] = geom[c][i] - geom[b][i]; 
    }
    for(int i = 0; i < DIM; i++) 
        v3[i] = cross(i, v1, v2);
    ans = dot(v0, v3)/sin(angle(a, b, c))/bond(a, b)/bond(c, b)/bond(d, b);
    if(ans < -1.0) ans = asin(-1.0);
    else if(ans > 1.0) ans = asin(1.0);
    else ans = asin(ans);
    return ans;
}

double Molecule::torsion(int a, int b, int c, int d)
{
    double v1[DIM], v2[DIM], v3[DIM], ans, mid1[DIM], mid2[DIM], sig[DIM];
    for(int i = 0; i < DIM; i++)
    {
        v1[i] = geom[b][i] - geom[a][i];
        v2[i] = geom[c][i] - geom[b][i]; 
        v3[i] = geom[d][c] - geom[c][i];
    }
    for(int i = 0; i < DIM; i++) 
    {
        mid1[i] = cross(i, v1, v2);
        mid2[i] = cross(i, v2, v3);
    }
    ans = dot(mid1, mid2)/bond(a, b)/bond(b, c)/bond(b, c)/bond(c, d)/sin(angle(a, b, c))/sin(angle(b, c, d));
    if(ans < -1.0) ans = acos(-1.0);
    else if(ans > 1.0) ans = acos(1.0);
    else ans = acos(ans);
    for(int i = 0; i < DIM; i++)
        sig[i] = cross(i, mid1, mid2);
    if(dot(sig, v2) < 0.0) ans = 0.0 - ans; 
    return ans;
}

double Molecule::com(int i)
{
    double ans = 0.0, mmass = 0.0;
    for(int l = 0; l < natom; l++)
    {
        ans += masses[zvals[l]] * geom[l][i];
        mmass += masses[zvals[l]];
    }
    ans = ans / mmass;
    return ans;
}

//Project 2
void Molecule::read_mw_hessian(const char *filename)
{
    FILE *fhessian = fopen(filename, "r");
    H.resize(3 * natom, 3 * natom);
    double sqrt_mass = 1.0;
    for(int i = 0; i < natom; i++)
        for(int j = 0; j < 3; j++)
            for(int _i = 0; _i < natom; _i++)
            {
                fscanf(fhessian, "%lf %lf %lf", &H(3*i + j, 3*_i) , &H(3*i + j, 3*_i + 1) , &H(3*i + j, 3*_i + 2));
                sqrt_mass = sqrt(masses[zvals[i]] * masses[zvals[_i]]);                   
                for(int _j = 0; _j < 3; _j++)
                {
                    H(3*i + j,3* _i + _j) /= sqrt_mass; 
                }
            }
    fclose(fhessian);

    return;
}




 
Molecule::Molecule(const char *filename)
{
    // open filename
    FILE *fmol = fopen(filename, "r");

    // read the number of atoms and the charge from filename
    fscanf(fmol, "%d %d", &natom, &charge);

    // allocate space
    zvals = new int[natom];
    geom = new double* [natom];
    for(int i = 0; i < natom; i++)
        geom[i] = new double[DIM];

    for(int i=0; i < natom; i++)
        fscanf(fmol, "%d %lf %lf %lf", zvals + i, geom[i] ,geom[i] + 1 ,geom[i] + 2);

    fclose(fmol);
}

Molecule::~Molecule()
{
    delete[] zvals;
    for(int i = 0; i < natom; i++)
        delete[] geom[i];
    delete[] geom;
}

double dot(double *v1, double *v2)
{
    double ans = 0;
    for(int i = 0; i < DIM; i++)
        ans += v1[i] * v2[i];
    return ans;
}

double cross(int i, double *v1, double *v2)
{
    double ans = 0;
    switch(i)
    {
        case 0:
            ans = v1[1] * v2[2] - v1[2] * v2[1];
            break;
        case 1:
            ans = v1[2] * v2[0] - v1[0] * v2[2];
            break;
        case 2:
            ans = v1[0] * v2[1] - v1[1] * v2[0];
            break;
    }
    return ans;
}