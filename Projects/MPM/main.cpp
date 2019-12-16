#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <random>
#include <chrono>

#include "SimulationDriver.h"
#include "PointCreation.h"
#include "TetToPoly.h"

#define CUSTOM_MESH 0

int main(int argc, char* argv[])
{
    using T = float;
    constexpr int dim = 3;
    using TV = Eigen::Matrix<T, dim, 1>;
    using TM = Eigen::Matrix<T,dim,dim>;
    using iV = Eigen::Matrix<int,dim,1>;

    // setup grid parameters
    // don't change min corner, it's assumed in weight computations
    TV gridMin = TV(0, 0, 0);
    TV gridMax = TV(1, 1, 1);
    T gridDx = 0.02;
    iV gridRes;
    gridRes[0] = ((gridMax[0] - gridMin[0]) / gridDx) + 1;
    gridRes[1] = ((gridMax[1] - gridMin[1]) / gridDx) + 1;
    gridRes[2] = ((gridMax[2] - gridMin[2]) / gridDx) + 1;

    int numGridValues = gridRes[0] * gridRes[1] * gridRes[2];

    // setup particle attributes
    T E = 10000;
    T nu = 0.3;
    T mu = E / (2.0 * (1.0 + nu));
    T lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    T rho = 10000;

    // sample particles
    std::vector<TV> startingPoints = PointCreation<T, dim>::createPoints(2, gridRes, gridMin, gridMax, false);

#if CUSTOM_MESH
    
#else // #if CUSTOM_MESH
    TV boxMin = TV(0.3, 0.1, 0.1);
    TV boxMax = TV(0.5, 0.3, 0.3);
    std::vector<TV> xp = PointCreation<T, dim>::selectInBox(startingPoints, boxMin, boxMax);
#endif // #else // #if CUSTOM_MESH

    int numParticles = xp.size();
    std::cout << "numParticles " << numParticles << std::endl;

    // particle initial volume
    std::vector<T> Vp0;
    Vp0.resize(numParticles, gridDx * gridDx / 4.0);

    // particle mass
    std::vector<T> mp;
    mp.resize(numParticles, Vp0[0] * rho); // fine bc all volumes are same

    // particle veloctiy
    std::vector<TV> vp;
    vp.resize(numParticles, TV::Zero());

    // particle deformation gradient
    std::vector<TM> Fp;
    Fp.resize(numParticles, TM::Identity());

    SimulationDriver<T,dim> driver;

    // simulate
    driver.mpm.numParticles = numParticles;
    driver.mpm.mp = mp;
    driver.mpm.vp = vp;
    driver.mpm.xp = xp;
    driver.mpm.Vp0 = Vp0;
    driver.mpm.Fp = Fp;
    driver.mpm.numGridValues = numGridValues;
    driver.mpm.gridRes = gridRes;
    driver.mpm.gridDx = gridDx;
    driver.mpm.mu = mu;
    driver.mpm.lambda = lambda;

    driver.run(480);

    return 0;
}
