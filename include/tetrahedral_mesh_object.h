#ifndef __TETRAHEDRAL_MESH_OBJECT__
#define __TETRAHEDRAL_MESH_OBJECT__

#include "mesh_object.h"

class TetrahedralMeshObject: public MeshObject {
    private:
        void set_mesh();
        void set_anchor_points();
        void set_springs();
        void set_massmatrix();

        double k;
        double size;
        char const* file;
    public:
        TetrahedralMeshObject(char const* file, double k=5.0, double size=1.0);

        Eigen::MatrixXi T;
};

#endif //__TETRAHEDRAL_MESH_OBJECT__