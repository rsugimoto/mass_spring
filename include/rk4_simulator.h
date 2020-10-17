#ifndef __RK4_SIMULATOR__
#define __RK4_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include "simulator.h"

//RK4 Integration Simulator
class RK4Simulator: public Simulator{
    private:
        Eigen::SparseMatrix<double> M_inverse;
    public:
        RK4Simulator(const MeshObject& obj, double dt);
        void forward_one_step();
};

#endif //__RK4_SIMULATOR__