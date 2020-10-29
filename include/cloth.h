#ifndef __CLOTH__
#define __CLOTH__

#include "mesh_object.h"

class Cloth: public MeshObject {
    private:
        void set_mesh();
        void set_anchor_points();
        void set_springs();

        double scale;
        int decimate;
    public:
        Cloth(double scale=1.0, int decimate=0);
};

#endif //__CLOTH__