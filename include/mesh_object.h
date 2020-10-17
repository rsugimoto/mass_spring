#ifndef __MESH_OBJECT__
#define __MESH_OBJECT__

#include <Eigen/Core>

class MeshObject {
    protected:
        virtual void set_mesh() = 0;
        virtual void set_anchor_points() = 0;
        virtual void set_springs() = 0;
        virtual void set_mass();
        void init();

    public:
        virtual ~MeshObject();
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd M; //mass
        Eigen::VectorXi A; //Anchor/fixed points
        Eigen::MatrixXi S; //Springs
        double k; //Spring Constant
};

#endif //__MESH_OBJECT__