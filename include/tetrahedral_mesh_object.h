#ifndef __TETRAHEDRAL_MESH_OBJECT__
#define __TETRAHEDRAL_MESH_OBJECT__

#include "mesh_object.h"

class TetrahedralMeshObject: public MeshObject {
    private:
        void set_mesh();
        void set_anchor_points();
        void set_springs();
        void set_mass();

        double scale;
        char const* file;
        int decimate;
    public:
        TetrahedralMeshObject(char const* file, double scale=1.0, int decimate=0);

        Eigen::MatrixXi T;
};

#endif //__TETRAHEDRAL_MESH_OBJECT__