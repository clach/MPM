#include <Eigen/Sparse>

#include <sys/stat.h>
#include <iostream>
#include "MPM.h"

#define SPHERE_COLLISION 1
#define FLOOR_COLLISION 1
#define DEBUG_MOMENTUM 1

template<class T, int dim>
class SimulationDriver {
public:
    using TV = Eigen::Matrix<T,dim,1>;
    using TM = Eigen::Matrix<T,dim,dim>;
    using iV = Eigen::Matrix<int,dim,1>;
    using SpMat = Eigen::SparseMatrix<T>;
    using Vec = Eigen::Matrix<T,Eigen::Dynamic,1>;

    MPM<T,dim> mpm;
    T dt;
    TV gravity;

    SimulationDriver()
    {
        gravity.setZero();
        gravity(1) = -9.8;
        dt = 0.0015;
    }

    void run(const int max_frame)
    {
         for (int frame = 1; frame < max_frame; frame++) {
            std::cout << "Frame " << frame << std::endl;

            int N_substeps = (int)(((T)1 / 24) / dt);
            for (int step = 1; step <= N_substeps; step++) {
                std::cout << "Step " << step << std::endl;
                
                // TODO: how often to run sim?

                // clear grid data
                mpm.clearGridData();

#if DEBUG_MOMENTUM
                TV Lp = mpm.computeParticleMomentum();
                std::cout << "Part momentum before P2G " << Lp(0) << " " << Lp(1) << " " << Lp(2) << std::endl;
#endif
                // P2G 
                // Splat mass and momentum, then find velocity for each grid node
                std::vector<int> activeNodes;
                mpm.transferP2G(mpm.xp, mpm.mp, mpm.vp, mpm.mg, mpm.vgn, activeNodes);
                //std::cout << "activeNodes.size " << activeNodes.size() << std::endl;

#if DEBUG_MOMENTUM
                TV Lg = mpm.computeGridMomentum(mpm.mg, mpm.vgn);
                std::cout << "Grid momentum after P2G " << Lg(0) << " " << Lg(1) << " " << Lg(2) << std::endl;
#endif
                // compute force
                mpm.addGravity(mpm.force, mpm.mg, activeNodes, gravity); 
                mpm.addElasticity(mpm.force, mpm.xp, mpm.Fp, mpm.Vp0); 

                // update velocity
                mpm.updateGridVelocity(mpm.mg, mpm.vgn, mpm.force, activeNodes, dt, mpm.vg);
            
                // boundary conditions
                mpm.setBoundaryVelocity(3, mpm.vg);      

#if DEBUG_MOMENTUM
                Lg = mpm.computeGridMomentum(mpm.mg, mpm.vg);
                std::cout << "Grid momentum before G2P " << Lg(0) << " " << Lg(1) << " " << Lg(2) << std::endl;
#endif                
                mpm.evolveF(dt, mpm.vg, mpm.xp, mpm.Fp);

                // G2P (including particle advection)
                mpm.transferG2P(dt, mpm.vgn, mpm.vg, 0.95, mpm.xp, mpm.vp);

#if DEBUG_MOMENTUM
                Lp = mpm.computeParticleMomentum();
                std::cout << "Part momentum after G2P " << Lp(0) << " " << Lp(1) << " " << Lp(2) << std::endl;
#endif                
            }

            mkdir("output/", 0777);
            std::string filename = "output/mpm" + std::to_string(frame) + ".poly";
            mpm.dumpPoly(filename);
            std::cout << std::endl;
        }
    }

};
