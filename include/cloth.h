#ifndef __CLOTH__
#define __CLOTH__

#include "mesh_object.h"

class Cloth: public MeshObject {
    private:
        void set_mesh();
        void set_anchor_points();
        void set_springs();

        double k;
        double size;
    public:
        Cloth(double k=5.0, double size=1.0);
};

#endif //__CLOTH__