#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

#include <vector>
#include <string>
#include <fstream>

#include "GridInterpolation.h"

template<class T, int dim>
class MPM {
public:
    using TV = Eigen::Matrix<T, dim, 1>;
    using TM = Eigen::Matrix<T, dim, dim>;
    using iV = Eigen::Matrix<int, dim, 1>;
    
    // particle data
    int numParticles;
    std::vector<TV> xp;
    std::vector<T> Vp0; 
    std::vector<T> mp;
    std::vector<TV> vp;
    std::vector<TM> Fp;

    // grid data
    int numGridValues;
    iV gridRes;
    T gridDx;
    std::vector<T> mg;
    std::vector<TV> vg;
    std::vector<TV> vgn;
    std::vector<TV> force;

    T mu;
    T lambda;

    MPM()
    {}

    void clearGridData() 
    {
        mg.clear();
        mg.resize(numGridValues, 0);
        vgn.clear();
        vgn.resize(numGridValues, TV::Zero());
        vg.clear();
        vg.resize(numGridValues, TV::Zero());
        force.clear();
        force.resize(numGridValues, TV::Zero());
    }

    TV computeParticleMomentum() const 
    {
        // compute total particle momentum
        TV pp = TV::Zero();
        for (int i = 0; i < numParticles; i++) 
        {
            pp += mp[i] * vp[i];
        }
        return pp;
    }

    TV computeGridMomentum() const
    {
        // compute total grid momentum
        TV pg = TV::Zero();
        for (int i = 0; i < numGridValues; i++) 
        {
            pg += mg[i] * vg[i];
        }
        return pg;
    }
    
    void transferP2G(const std::vector<TV>& xp, const std::vector<T>& mp, const std::vector<TV>& vp, 
                    std::vector<T>& mg, std::vector<TV>& vgn, std::vector<int>& activeNodes) 
    {
        for (int p = 0; p < numParticles; p++)
        {
            TV X = xp[p];
            iV index;
            index[0] = (int)(X[0] / gridDx);
            index[1] = (int)(X[1] / gridDx);
            index[2] = (int)(X[2] / gridDx);

            TV wX;
            TV wY;
            TV wZ;
            int baseNodeX = GridInterpolation<T, dim>::computeWeights1D(index[0], wX);
            int baseNodeY = GridInterpolation<T, dim>::computeWeights1D(index[1], wY);
            int baseNodeZ = GridInterpolation<T, dim>::computeWeights1D(index[2], wZ);

            for (int i = 0; i < dim; i++) 
            {
                T w_i = wX[i];
                int node_i = baseNodeX + i;

                for (int j = 0; j < dim; j++) 
                {
                    T w_ij = w_i * wY[j];
                    int node_j = baseNodeY + j;

                    for (int k = 0; k < dim; k++) 
                    {
                        T w_ijk = w_ij * wZ[k];
                        int node_k = baseNodeZ + k;

                        // splat mass
                        int gridIndex = node_i + node_j * gridRes[0] + node_k * gridRes[0] * gridRes[1];
                        mg[gridIndex] += mp[p] * w_ijk;

                        // splat momentum (store in vg temporarily)
                        vgn[gridIndex] += (w_ijk * mp[p]) * vp[p];
                    }
                }
            }
        }

        // determine which nodes are active (have nearby particles i.e. nonzero mass)
        for (int i = 0; i < numGridValues; i++) 
        {
            if (mg[i] != 0) // don't need to do FP comparison
            {
                activeNodes.push_back(i);
                vgn[i] = vg[i] / mg[i]; // remember we stored momentum in vg 
            } 
            else 
            {
                vgn[i] = TV::Zero();
            }
        }
    }

    void transferG2P(T dt, const std::vector<TV>& vgn, const std::vector<TV>& vg, T flip,
                    std::vector<TV>& xp, std::vector<TV>& vp) 
    {
        for (int p = 0; p < numParticles; p++)
        {
            TV X = xp[p];
            iV index;
            index[0] = (int)(X[0] / gridDx);
            index[1] = (int)(X[1] / gridDx);
            index[2] = (int)(X[2] / gridDx);

            TV wX, wY, wZ;
            int baseNodeX = GridInterpolation<T, dim>::computeWeights1D(index[0], wX);
            int baseNodeY = GridInterpolation<T, dim>::computeWeights1D(index[1], wY);
            int baseNodeZ = GridInterpolation<T, dim>::computeWeights1D(index[2], wZ);

            // compute v_pic and v_flip
            TV v_pic = TV::Zero();
            TV v_flip = vp[p];

            for (int i = 0; i < dim; i++) 
            {
                T w_i = wX[i];
                int node_i = baseNodeX + i;

                for (int j = 0; j < dim; j++) 
                {
                    T w_ij = w_i * wY[j];
                    int node_j = baseNodeY + j;

                    for (int k = 0; k < dim; k++) 
                    {
                        T w_ijk = w_ij * wZ[k];
                        int node_k = baseNodeZ + k;

                        int gridIndex = node_i + node_j * gridRes[0] + node_k * gridRes[0] * gridRes[1];
                        v_pic += w_ijk * vg[gridIndex];
                        v_flip += w_ijk * (vg[gridIndex] - vgn[gridIndex]);
                    }
                }
            }

            // update particle velocity and particle position
            vp[p] = (1.0 - flip) * v_pic + flip * v_flip;
            xp[p] += dt * v_pic; // only use pic when updating position
        }
    }

    void addGravity(std::vector<TV>& force, const std::vector<T>& mg, 
                    const std::vector<int>& activeNodes, const TV& gravity) 
    {
        // add gravity to force
        for (size_t i = 0; i < activeNodes.size(); i++) {
            force[activeNodes[i]] += mg[activeNodes[i]] * gravity;
        }
    }

    void polarSVD(const TM& F, TM& U, TM& Sigma, TM& V) 
    {
        // compute regular SVD
        Eigen::JacobiSVD<TM> svd(F, Eigen::ComputeThinV | Eigen::ComputeThinU);
        U = svd.matrixU();
        Sigma = svd.singularValues().asDiagonal();
        V = svd.matrixV();

        // convert to Polar SVD
        if (U.determinant() < 0) 
        {
            U(0, 2) *= -1;
            U(1, 2) *= -1;
            U(2, 2) *= -1;
            Sigma(2, 2) *= -1;
        }

        if (V.determinant() < 0) 
        {
            V(0, 2) *= -1;
            V(1, 2) *= -1;
            V(2, 2) *= -1;
            Sigma(2, 2) *= -1;
        }
    }

    TM fixCorotated(const TM& F) 
    {
        TM U;
        TM Sigma;
        TM V;
        polarSVD(F, U, Sigma, V);
        TM R = U * V.transpose();
        TM newF = U * Sigma * V.transpose();

        float J = newF.determinant();

        // A = (JF)^-T
        TM A = newF.adjoint().transpose();

        // compute P = dPsi / dF
        TM P = 2.f * mu * (newF - R) + lambda * (J - 1.f) * A;

        return P;
    }

    void addElasticity(std::vector<TV>& force, const std::vector<TV>& xp, 
                    const std::vector<TM>& Fp, const std::vector<T>& Vp0) 
    {
        for (int p = 0; p < numParticles; p++) 
        {
            TM currFp = Fp[p];
            TM P = fixCorotated(currFp);
            TM Vp0_P_FTrans = Vp0[p] * P * currFp.transpose();

            TV X = xp[p];
            iV index = iV((int)(X(0) / gridDx), (int)(X(1) / gridDx), (int)(X(2) / gridDx));

            TV wX, wY, wZ, dwX, dwY, dwZ;
            int baseNodeX = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[0], wX, dwX);
            int baseNodeY = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[1], wY, dwY);
            int baseNodeZ = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[2], wZ, dwZ);

            for (int i = 0; i < dim; i++) 
            {
                T w_i = wX(i);
                T dw_i_dx_i = dwX(i) / gridDx;
                int node_i = baseNodeX + i;

                for (int j = 0; j < dim; j++) 
                {
                    T w_j = wY(j);
                    T w_ij = w_i * w_j;
                    T dw_ij_dx_i = dw_i_dx_i * w_j;
                    T dw_ij_dx_j = dwY(j) / gridDx * w_i;

                    int node_j = baseNodeY + j;

                    for (int k = 0; k < dim; k++) 
                    {
                        T w_k = wZ(k);
                        T dw_ijk_dx_i = dw_ij_dx_i * w_k;
                        T dw_ijk_dx_j = dw_ij_dx_j * w_k;
                        T dw_ijk_dx_k = dwZ(k) / gridDx * w_ij;

                        int node_k = baseNodeZ + k;

                        TV grad_w = TV(dw_ijk_dx_i, dw_ijk_dx_j, dw_ijk_dx_k);   

                        int gridIndex = node_i + node_j * gridRes[0] + node_k * gridRes[0] * gridRes[1];
                        force[gridIndex] += Vp0_P_FTrans * grad_w;
                    }
                }
            }
        }
    }

    void updateGridVelocity(const std::vector<T>& mg, const std::vector<TV>& vgn, 
                            const std::vector<TV>& force, const std::vector<int>& activeNodes,
                            T dt, std::vector<TV>& vg) 
    {
        // update velocity from force
        for (size_t i = 0; i < activeNodes.size(); i++) 
        {   
            vg[activeNodes[i]] = vgn[activeNodes[i]] + dt * force[activeNodes[i]] / mg[activeNodes[i]];
        }
    }

    void setBoundaryVelocity(int thickness, std::vector<TV>& vg)
    {
        // set domain boundary velocities

        // TODO
        // min x direction
        for (int i = 0; i < thickness; i++) 
        {
            for (int j = 0; j < gridRes[1]; j++)
            {
                for (int k = 0; k < gridRes[2]; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }

        // max x direction
        for (int i = gridRes[0] - thickness; i < gridRes[0]; i++) 
        {
            for (int j = 0; j < gridRes[1]; j++)
            {
                for (int k = 0; k < gridRes[2]; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }

        // min y direction
        for (int i = 0; i < gridRes[0]; i++) 
        {
            for (int j = 0; j < thickness; j++)
            {
                for (int k = 0; k < gridRes[2]; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }

        // max y direction
        for (int i = 0; i < gridRes[0]; i++) 
        {
            for (int j = gridRes[1] - thickness; j < gridRes[1]; j++)
            {
                for (int k = 0; k < gridRes[2]; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }

        // min z direction
        for (int i = 0; i < gridRes[0]; i++) 
        {
            for (int j = 0; j < gridRes[1]; j++)
            {
                for (int k = 0; k < thickness; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }

        // max z direction
        for (int i = 0; i < gridRes[0]; i++) 
        {
            for (int j = 0; j < gridRes[1]; j++)
            {
                for (int k = gridRes[2] - thickness; k < gridRes[2]; k++)
                {
                    int index = i + j * gridRes[0] + k * gridRes[0] * gridRes[1];
                    vg[index] = TV::Zero();
                }
            }
        }
    }

    void evolveF(T dt, const std::vector<TV>& vg, const std::vector<TV>& xp, std::vector<TM>& Fp ) 
    {
        for (int p = 0; p < numParticles; p++) 
        {
            TM currFp = Fp[p];

            // compute grad vp
            TM grad_vp = TM::Zero();

            TV X = xp[p];
            iV index = iV((int)(X[0] / gridDx), (int)(X[1] / gridDx), (int)(X[2] / gridDx));

            TV wX, wY, wZ, dwX, dwY, dwZ;
            int baseNodeX = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[0], wX, dwX);
            int baseNodeY = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[1], wY, dwY);
            int baseNodeZ = GridInterpolation<T, dim>::computeWeightsWithGradient1D(index[2], wZ, dwZ);

            for (int i = 0; i < dim; i++) 
            {
                T w_i = wX(i);
                T dw_i_dx_i = dwX(i) / gridDx;
                int node_i = baseNodeX + i;

                for (int j = 0; j < dim; j++) 
                {
                    T w_j = wY(j);
                    T w_ij = w_i * w_j;
                    T dw_ij_dx_i = dw_i_dx_i * w_j;
                    T dw_ij_dx_j = dwY(j) / gridDx * w_i;

                    int node_j = baseNodeY + j;

                    for (int k = 0; k < dim; k++) 
                    {
                        T w_k = wZ(k);
                        T dw_ijk_dx_i = dw_ij_dx_i * w_k;
                        T dw_ijk_dx_j = dw_ij_dx_j * w_k;
                        T dw_ijk_dx_k = dwZ(k) / gridDx * w_ij ;

                        int node_k = baseNodeZ + k;

                        TV grad_w = TV(dw_ijk_dx_i, dw_ijk_dx_j, dw_ijk_dx_k);   

                        int gridIndex = node_i + node_j * gridRes[0] + node_k * gridRes[0] * gridRes[1];
                        
                        grad_vp += vg[gridIndex] * grad_w.transpose();          
                    }
                }
            }

            Fp[p] = (TM::Identity() + dt * grad_vp) * currFp;     
        }
    }

    void dumpPoly(std::string filename)
    {
        std::ofstream fs;
        fs.open(filename);
        fs << "POINTS\n";
        int count = 0;
        for (auto X : xp) {
            fs << ++count << ":";
            for (int i = 0; i < dim; i++)
                fs << " " << X(i);
            if (dim == 2)
                fs << " 0";
            fs << "\n";
        }
        fs << "POLYS\n";
        count = 0;
        //for (const Eigen::Matrix<int, 2, 1>& seg : segments)
            //fs << ++count << ": " << seg(0) + 1 << " " << seg(1) + 1 << "\n"; // poly segment mesh is 1-indexed
        fs << "END\n";
        fs.close();
    }
    
};
