#ifndef __DTTLE_SIMULATOR__
#define __DTTLE_SIMULATOR__

#include "simulator.h"

//Discrete-time Toda Latice Equation Simulator
class DTTLESimulator: public Simulator{
    private:
        const double a;
        const double b;
        
    public:
        DTTLESimulator(const MeshObject& obj, double a, double b, double dt);
        void step();
};

#endif //__DTTLE_SIMULATOR__