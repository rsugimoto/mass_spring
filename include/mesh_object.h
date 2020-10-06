#ifndef __MESH_OBJECT__
#define __MESH_OBJECT__

#include <Eigen/Core>

class MeshObject {
    protected:
        virtual void set_mesh() = 0;
        virtual void set_anchor_points() = 0;
        virtual void set_springs() = 0;
        void set_rest_lengths();
        virtual void set_massmatrix();
        void init();

    public:
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        Eigen::VectorXd RL; //rest length
        Eigen::VectorXd M; //mass matrix
        Eigen::VectorXi A;
        Eigen::MatrixXi S; //Springs
        Eigen::VectorXd K; //Spring Constants
};

#endif //__MESH_OBJECT__