#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <fstream>

template<class T, int dim>
class GridInterpolation {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using iV = Eigen::Matrix<int,dim,1>;

    static int computeWeights1D(int x, TV& w) 
    {
        // compute 1D quadratic B-spline weights
        // x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
        // w is 1x3 (row vector)
        int baseNode = floor(x - 0.5) + 1;

        w = TV::Zero();

        int d0 = x - baseNode + 1;
        T z = 1.5 - (T)d0;
        T z2 = z * z;
        w[0] = 0.5 * z2;

        int d1 = d0 - 1;
        w[1] = 0.75 - d1 * d1;

        int d2 = 1 - d1;
        T zz = 1.5 - (T)d2;
        T zz2 = zz * zz;
        w[2] = 0.5 * zz2;

        return baseNode;
    }

    static int computeWeightsWithGradient1D(int x, TV& w, TV& dw) 
    {
        // compute 1D quadratic B-spline weights
        // x is assumed to be scaled in the index space (i.e., it is in a dx=1 grid)
        // w is 1x3 (row vector)
        int baseNode = floor(x - 0.5) + 1;

        w = TV::Zero();
        dw = TV::Zero();

        int d0 = x - baseNode + 1;
        T z = 1.5 - (T)d0;
        T z2 = z * z;
        w[0] = 0.5 * z2;

        int d1 = d0 - 1;
        w[1] = 0.75 - d1 * d1;

        int d2 = 1 - d1;
        T zz = 1.5 - (T)d2;
        T zz2 = zz * zz;
        w[2] = 0.5 * zz2;

        dw(1) = -z;
        dw(2) = -2.0 * d1;
        dw(3) = zz;

        return baseNode;
    }
};