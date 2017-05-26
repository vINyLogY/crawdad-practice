#include <iostream>
#include <string>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

#define PI 3.1415926538979
#define DIM 3


using namespace std;
 
class Molecule
{
    public:
        int natom;
        int charge;
        int *zvals;
        double **geom;
        string point_group;
        Matrix H;
        
        //Project 1
        void print_geom();
        void rotate_z(double phi);
        void translate(double x, double y, double z);
        double bond(int atom1, int atom2);
        double angle(int atom1, int atom2, int atom3);
        double angle_op(int atom1, int atom2, int atom3, int atom4);
        double torsion(int atom1, int atom2, int atom3, int atom4);
        double com(int i);

        //Project 2
        void read_mw_hessian(const char *filename);        

        Molecule(const char *filename);
        ~Molecule();
};

double dot(double *v1, double *v2);
double cross(int i, double *v1, double *v2);
