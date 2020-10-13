#ifndef __LITI_SIMULATOR__
#define __LITI_SIMULATOR__

//Reference
//https://github.com/dilevin/CSC417-a2-mass-spring-3d

#include "simulator.h"

//Linearly-Implicit Time Integration Simulator
class LITISimulator: public Simulator{
    public:
        LITISimulator(const MeshObject& obj, double dt);
        void forward_one_step();
};

#endif //__LITI_SIMULATOR__