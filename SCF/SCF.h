#include <iostream>
#include <string>

#include "Eigen/Dense"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class Basis
{
    public:
        Matrix S, T, V, H, Sx;
        Matrix Fx, Eps, Cx, C, D;
        int max, occ;
        double Enuc, *tei, Etot, Eele;

        void set_fx(int i);
        void set_eps(int i);
        void set_d(int i);

        Basis(const char *enuc, const char *s, const char *t, const char *v, const char *eri);
        ~Basis();

        
};

Matrix fr_basis(const char *filename);






