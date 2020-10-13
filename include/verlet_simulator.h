#ifndef __VERLET_SIMULATOR__
#define __VERLET_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include "simulator.h"

//Linearly-Implicit Time Integration Simulator
class VerletSimulator: public Simulator{
    private:
        Eigen::SparseMatrix<double> M_inverse;
    public:
        void set_initial_config();
        VerletSimulator(const MeshObject& obj, double dt);
        void forward_one_step();
};

#endif //__VERLET_SIMULATOR__