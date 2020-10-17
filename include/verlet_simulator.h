#ifndef __VERLET_SIMULATOR__
#define __VERLET_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include "simulator.h"

//Verlet Integration Simulator
class VerletSimulator: public Simulator{
    private:
        Eigen::SparseMatrix<double> M_inverse;
    public:
        VerletSimulator(const MeshObject& obj, double dt);
        void set_initial_config();
        void forward_one_step();
};

#endif //__VERLET_SIMULATOR__