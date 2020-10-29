#ifndef __RK4_SIMULATOR__
#define __RK4_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include "simulator.h"

//RK4 Integration Simulator
class RK4Simulator: public Simulator{
    public:
        RK4Simulator(const MeshObject& obj, double k, double dt);
        void step();
};

#endif //__RK4_SIMULATOR__