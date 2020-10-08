#ifndef __LTITI_SIMULATOR__
#define __LTITI_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <mesh_object.h>

//Linearly Time Implicit Time Integration Simulator
class LTITISimulator{
    public:
        LTITISimulator(const MeshObject& obj, double dt);
        void forward_one_step();
        Eigen::Map<const Eigen::MatrixXd> V(); //behaves similarly to a constant reference
        void set_initial_config();
    private:
        Eigen::SparseMatrix<double> M;
        Eigen::VectorXd q_curr, q_prev, q_dot_curr, q_dot_prev, L;
        Eigen::SparseMatrix<double> P;

        const MeshObject& obj;
        void d2Vdq2(Eigen::SparseMatrix<double>& d2Vdq2);
        void dVdq(Eigen::VectorXd& dVdq);
        double dt;
};

#endif //__LTITI_SIMULATOR__