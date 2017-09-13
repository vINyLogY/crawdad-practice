#include <iostream>
#include <string>

#include "Eigen/Dense"

#define INDEX(i, j) ((i > j) ? (((i)*(i+1)/2)+j) : ((j)*(j+1)/2+i))

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class Basis
{
    public:
        Matrix S, T, V, H, Sx;
        Matrix F, Fx, Eps, Cx, C, D;
        Matrix _D; //储存上一次迭代的内容
        int max, occ, iter;
        double Enuc, *tei, Etot, Eele;
        double _Eele;

        void set_fx();
        void set_eps();
        void set_d();

        bool e_conv(double delta);
        bool rms_conv(double delta);

        void iteration(int oao, double err);

        Matrix MOF();
        double Mu(const char *mu);

        Basis(const char *enuc, const char *s, const char *t, const char *v, const char *eri);
        ~Basis();   
};

Matrix fr_basis(const char *filename);








